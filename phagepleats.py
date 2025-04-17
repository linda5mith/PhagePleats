import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
import faiss
from sklearn.metrics import jaccard_score
from collections import defaultdict
from joblib import Parallel, delayed
import time

# 1. Read Foldseek result file
def read_search(path_to_foldseek_search):
    """
    Reads Foldseek search results and removes '.pdb' extensions from query and target columns.

    Args:
        path_to_foldseek_search (str): Path to the Foldseek result file (TSV).

    Returns:
        pd.DataFrame: Cleaned Foldseek search DataFrame with column names.
    """
    search = pd.read_csv(path_to_foldseek_search, sep='\t', header=None)
    search.columns = ['query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits']
    search['query'] = search['query'].str.replace('.pdb', '', regex=True)
    search['target'] = search['target'].str.replace('.pdb', '', regex=True)
    return search

# 2. Read metadata mapping protein -> genome
def read_metadata(path_to_input_metadata):
    """
    Reads metadata that maps proteins to genomes.

    Args:
        path_to_input_metadata (str): Path to the metadata CSV.

    Returns:
        pd.DataFrame: Metadata DataFrame with 'genome' and 'protein' columns.
    """
    return pd.read_csv(path_to_input_metadata)

# 3. Map each query to its genome using metadata
def map_query_to_genome(search, metadata):
    """
    Adds a 'query_genome' column to the search DataFrame by mapping each query to its genome using metadata.

    Args:
        search (pd.DataFrame): Foldseek result DataFrame.
        metadata (pd.DataFrame): DataFrame with columns 'genome' and 'protein'.

    Returns:
        pd.DataFrame: Updated search DataFrame with 'query_genome' column.
    """
    protein_to_genome = dict(zip(metadata['protein'], metadata['genome']))
    search['query_genome'] = search['query'].map(protein_to_genome)
    return search

def create_input_matrix(search, presence_absence_path):
    """Creates a presence/absence matrix for input genomes based on detected clusters.

    Args:
        search (pd.DataFrame): Search DataFrame with 'query_genome' and 'target'.
        presence_absence_path (str): Path to presence/absence file with cleaned cluster names.

    Returns:
        pd.DataFrame: Presence/absence matrix (genomes x clusters).
    """
    print("Building presence/absence matrix...")
    presence_absence = pd.read_csv(presence_absence_path)
    presence_absence[['accession', 'function']] = presence_absence['cluster_ID_function'].str.split(":", n=1, expand=True)
    presence_absence.set_index('cluster_ID_function', inplace=True)

    # Make mapping from accession -> cluster_ID_function
    accession_to_cluster = presence_absence['accession'].to_dict()
    cluster_lookup = {v: k for k, v in accession_to_cluster.items()}

    # Filter only known clusters in the lookup
    search['matched_cluster'] = search['target'].map(cluster_lookup)
    valid = search.dropna(subset=['query_genome', 'matched_cluster'])

    # Pivot table into presence/absence matrix
    pivot = (
        valid
        .drop_duplicates(subset=['query_genome', 'matched_cluster'])
        .assign(present=1)
        .pivot(index='query_genome', columns='matched_cluster', values='present')
        .fillna(0)
        .astype('float32')
    )

    # Reindex to ensure all clusters in same order as original
    final = pivot.reindex(columns=presence_absence.index, fill_value=0)
    return final

# 5. Read in models and run prediction
def load_models_from_folder(folder_path):
    """
    Loads pickled models from a folder.

    Args:
        folder_path (str): Path to folder containing .pkl model files.

    Returns:
        dict: Dictionary of models keyed by taxonomic label.
    """
    models = {}
    for fname in os.listdir(folder_path):
        if fname.endswith(".pkl"):
            taxon = fname.replace(".pkl", "")
            full_path = os.path.join(folder_path, fname)
            with open(full_path, "rb") as f:
                models[taxon] = pickle.load(f)
    return models

def predict_from_models(X, models, rank):
    """
    Predicts taxonomic labels and probabilities using models.

    Args:
        X (pd.DataFrame): Input presence/absence matrix.
        models (dict): Dictionary of models.
        rank (str): Taxonomic rank.

    Returns:
        pd.DataFrame: Predicted labels and probabilities.
        pd.DataFrame: Full probability matrix.
    """
    prob_list = []

    for taxon, model in tqdm(models.items(), desc=f"Predicting for {rank}"):
        probs = model.predict_proba(X)[:, 1]  # probability of class 1
        prob_list.append(probs)

    prob_matrix = pd.DataFrame(np.column_stack(prob_list), index=X.index, columns=models.keys())
    pred_labels = prob_matrix.idxmax(axis=1)
    pred_scores = prob_matrix.max(axis=1)

    df_out = pd.DataFrame({
        f"{rank}": pred_labels,
        f"{rank}_prob": pred_scores
    }, index=X.index)

    df_out.index.name = 'Genome'
    df_out = df_out.sort_index()

    return df_out, prob_matrix

def predict_all_ranks(X, ranks=["Order", "Family", "Subfamily", "Genus"], model_base="data/models", out_dir="outputs"):
    """
    Predicts taxonomic ranks using saved models for each rank and saves results to the specified output directory.

    Args:
        X (pd.DataFrame): Input presence/absence matrix.
        ranks (list): List of taxonomic ranks.
        model_base (str): Base directory containing rank-specific model folders.
        out_dir (str): Output directory where results will be saved.

    Returns:
        pd.DataFrame: Final DataFrame containing predictions for all ranks.
    """
    all_preds = []
    os.makedirs(out_dir, exist_ok=True)
    prob_dir = os.path.join(out_dir, "rank_probabilities")
    os.makedirs(prob_dir, exist_ok=True)

    for rank in ranks:
        print(f"\n🔮🧬✨ Predicting for rank: {rank}")
        model_path = os.path.join(model_base, rank)
        if not os.path.exists(model_path):
            print(f"❌ No model directory for {rank}: {model_path}")
            continue

        models = load_models_from_folder(model_path)
        if not models:
            print(f"⚠️ No models found in {model_path}")
            continue

        df_rank, df_proba = predict_from_models(X, models, rank)
        df_rank.index = X.index
        all_preds.append(df_rank)
        df_proba = df_proba.sort_index()
        df_proba.to_csv(os.path.join(prob_dir, f"{rank}_proba.csv"))

    final_df = pd.concat(all_preds, axis=1)
    final_df.index.name = "Genome"
    final_df = final_df.sort_index()

    print(f"\n✅ Final predictions saved to '{os.path.join(out_dir, 'taxonomy_predictions.csv')}'")
    print(f"🧬📂 Full probability matrices saved to '{prob_dir}/'")

    return final_df

# 7. Compute distances to training matrix using FAISS
def compute_closest_training_genomes(presence_absence_path, input_matrix, taxonomy_df):
    """
    Finds the closest training genome for each input genome and computes:
    - Euclidean distance
    - % shared proteins
    - Taxonomic ranks of closest genome (Order, Family, Subfamily, Genus)

    Args:
        presence_absence_path (str): Path to the training presence/absence matrix CSV.
        input_matrix (pd.DataFrame): Presence/absence matrix of input genomes.
        taxonomy_df (pd.DataFrame): DataFrame with columns ['Leaves', 'Order', 'Family', 'Subfamily', 'Genus'].

    Returns:
        pd.DataFrame: DataFrame with closest genome info, distance, % shared, and taxonomic ranks.
    """
    # Load training matrix
    presence_absence = pd.read_csv(presence_absence_path, index_col=0, compression='gzip')
    training_matrix = presence_absence.astype('float32').T
    input_matrix = input_matrix.astype('float32')

    # Build FAISS index
    index = faiss.IndexFlatL2(training_matrix.shape[1])
    index.add(training_matrix.values)

    # Search
    distances, indices = index.search(input_matrix.values, k=1)
    closest_genomes = training_matrix.index[indices.flatten()]
    sqrt_distances = np.sqrt(distances.flatten())

    # Set taxonomy index
    taxonomy_df = taxonomy_df.set_index("Leaves").astype(str)

    # Compute results
    results = []
    for genome, closest, dist in zip(input_matrix.index, closest_genomes, sqrt_distances):
        input_vec = input_matrix.loc[genome].astype(bool)
        closest_vec = training_matrix.loc[closest].astype(bool)
        shared = jaccard_score(input_vec, closest_vec, average='binary', zero_division=0)

        tax = taxonomy_df.loc[closest].to_dict() if closest in taxonomy_df.index else {}

        results.append({
            "input_genome": genome,
            "closest_training_genome": closest,
            "euclidean_distance": dist,
            "%_shared_proteins": shared,
            "Order": tax.get("Order", np.nan),
            "Family": tax.get("Family", np.nan),
            "Subfamily": tax.get("Subfamily", np.nan),
            "Genus": tax.get("Genus", np.nan)
        })

    return pd.DataFrame(results).round(3).sort_values(by='input_genome')

# --- Helper: Preprocess matrices and build clade members ---
def preprocess_matrices(training_matrix, input_matrix, taxonomy_df, preds):
    training_matrix_f32 = training_matrix.astype('float32')
    training_matrix_bool = training_matrix.astype(bool)
    input_matrix_f32 = input_matrix.astype('float32')
    input_matrix_bool = input_matrix.astype(bool)

    clade_members = defaultdict(dict)
    for rank in ["Order", "Family", "Subfamily", "Genus"]:
        taxonomy_df[rank] = taxonomy_df[rank].astype(str)
        preds[rank] = preds[rank].astype(str)
        for clade, clade_df in taxonomy_df.groupby(rank):
            matched = clade_df.index.intersection(training_matrix.index)
            if len(matched) > 0:
                clade_members[rank][clade] = matched

    return training_matrix_f32, training_matrix_bool, input_matrix_f32, input_matrix_bool, clade_members

# --- Helper: Build FAISS indexes ---
def build_faiss_indexes(clade_members, training_matrix_f32):
    faiss_indexes = defaultdict(dict)
    for rank in clade_members:
        for clade, genomes in clade_members[rank].items():
            vecs = training_matrix_f32.loc[genomes].values
            if vecs.shape[0] > 0:
                index = faiss.IndexFlatL2(vecs.shape[1])
                index.add(vecs)
                faiss_indexes[rank][clade] = index
    return faiss_indexes

# --- Helper: Process a single genome ---
def process_genome(genome, input_matrix_f32, input_matrix_bool, training_matrix_bool, clade_members, faiss_indexes, preds):
    row = {"Genome": genome}
    input_vec = input_matrix_f32.loc[genome].values
    input_vec_bool = input_matrix_bool.loc[genome].values

    for rank in ["Order", "Family", "Subfamily", "Genus"]:
        pred_clade = preds.loc[genome, rank]
        matched_genomes = clade_members[rank].get(pred_clade, [])

        if len(matched_genomes) == 0:
            row[f"%_shared_with_predicted_{rank}"] = np.nan
            row[f"eucl_dist_to_predicted_{rank}"] = np.nan
        else:
            clade_bools = training_matrix_bool.loc[matched_genomes].values
            jaccard_vals = np.mean([
                jaccard_score(input_vec_bool, clade_bools[i], average='binary', zero_division=0)
                for i in range(clade_bools.shape[0])
            ])
            row[f"%_shared_with_predicted_{rank}"] = jaccard_vals

            faiss_index = faiss_indexes[rank].get(pred_clade)
            if faiss_index:
                D, _ = faiss_index.search(input_vec.reshape(1, -1), k=faiss_index.ntotal)
                row[f"eucl_dist_to_predicted_{rank}"] = np.sqrt(np.mean(D))
            else:
                row[f"eucl_dist_to_predicted_{rank}"] = np.nan

    return row

# --- Helper: Classify novelty based on z-scores ---
def classify_novelty(z_s, z_d, rank):
    if pd.isna(z_s) or pd.isna(z_d):
        return "Unknown"
    for threshold, label in zip([3, 2, 1], ["order", "family", "subfamily", "genus"]):
        if z_s <= -threshold or z_d >= threshold:
            return f"Potential new {label if rank.lower() == label else rank.lower()}"
    return "Likely member"

# --- Additional: Compute z-scores and flags for all ranks ---
def compute_z_scores_and_flag_all_ranks(phage_df, intra_df):
    for rank in ["Order", "Family", "Subfamily", "Genus"]:
        z_shared_list = []
        z_euclid_list = []
        flag_list = []

        for _, row in phage_df.iterrows():
            predicted_clade = row[rank]
            clade_metrics = intra_df[(intra_df['rank'] == rank) & (intra_df['clade'] == predicted_clade)]

            if clade_metrics.empty:
                z_shared = z_euclid = None
                novelty_flag = "Clade not in intra table"
            else:
                mean_shared = clade_metrics["intra_avg_shared_proteins"].values[0]
                std_shared = clade_metrics["intra_std_shared_proteins"].values[0]
                mean_euclid = clade_metrics["intra_avg_euclidean"].values[0]
                std_euclid = clade_metrics["intra_std_euclidean"].values[0]

                try:
                    shared = row[f"%_shared_with_predicted_{rank}"]
                    euclid = row[f"eucl_dist_to_predicted_{rank}"]
                    if std_shared == 0 or std_euclid == 0:
                        z_shared = z_euclid = np.nan
                    else:
                        z_shared = (shared - mean_shared) / std_shared
                        z_euclid = (euclid - mean_euclid) / std_euclid
                except:
                    z_shared = z_euclid = np.nan

                novelty_flag = classify_novelty(z_shared, z_euclid, rank)

            z_shared_list.append(z_shared)
            z_euclid_list.append(z_euclid)
            flag_list.append(novelty_flag)
            
            print(predicted_clade)
            print(z_shared)
            print(z_euclid)
            print(novelty_flag)
            print('-----------------------')

        phage_df[f"{rank}_z_shared_proteins"] = z_shared_list
        phage_df[f"{rank}_z_euclidean_distance"] = z_euclid_list
        phage_df[f"{rank}_novelty_flag"] = flag_list

    return phage_df


# #Novelty-aware post-prediction QC layer 👾🧬✨
def compute_clade_novelty_summary(presence_absence_path, input_matrix, taxonomy_df, preds, intra_rank_relatedness):
    print("🔍 Loading training presence/absence matrix...")
    presence_absence = pd.read_csv(presence_absence_path, index_col=0, compression='gzip')
    training_matrix = presence_absence.astype('float32').T
    input_matrix = input_matrix.astype('float32')
    taxonomy_df = taxonomy_df.set_index("Leaves")

    print("🔗 Finding closest training genome (Euclidean distance) to input genome...")
    index = faiss.IndexFlatL2(training_matrix.shape[1])
    index.add(training_matrix.values)
    distances, indices = index.search(input_matrix.values, k=1)
    closest_df = pd.DataFrame({
        'Genome': input_matrix.index,
        'closest_training_genome': training_matrix.index[indices.flatten()],
        'euclidean_dist_to_closest_hit': np.sqrt(distances.flatten())
    }).set_index("Genome")

    print("🔬 Computing % shared proteins + FAISS distance to predicted clades...")
    start = time.perf_counter()
    training_matrix_f32, training_matrix_bool, input_matrix_f32, input_matrix_bool, clade_members = preprocess_matrices(
        training_matrix, input_matrix, taxonomy_df, preds
    )
    faiss_indexes = build_faiss_indexes(clade_members, training_matrix_f32)

    results = []
    for genome in tqdm(input_matrix.index, desc="💫 Calculating per-genome clade similarity"):
        row = process_genome(genome, input_matrix_f32, input_matrix_bool, training_matrix_bool, clade_members, faiss_indexes, preds)
        results.append(row)

    end = time.perf_counter()
    print(f"⏱️ Took {end - start:.2f} seconds for intra-clade similarity calculations")

    clade_df = pd.DataFrame(results).set_index("Genome")
    final_df = closest_df.join(clade_df)

    print("🧬 Assigning novelty scores...")
    for rank in ["Order", "Family", "Subfamily", "Genus"]:
        final_df[rank] = preds[rank]

    final_df = compute_z_scores_and_flag_all_ranks(final_df, intra_rank_relatedness)
    final_df.index.name = "Genome"
    final_df = final_df.round(2).sort_index()
    print("✅ Novelty summary complete.")
    return final_df

# Main Snakemake execution
if __name__ == "__main__":
    def print_phage_ascii():
        print(r'''
          ___
        /     \
       | o   o |
        \ ˇ ˇ /
         |||||
         |||||
          |||
       ___|||___
      /   |||   \
    PhagePleats is predicting... 🧬🚀✨
              
    Uncovering the virosphere one phage at a time... 🔬🦧👾
    ''')
        
    # ... call it before prediction starts lol:
    print_phage_ascii()

    search_result = snakemake.input.search_result
    metadata = snakemake.input.query_metadata
    presence_absence = snakemake.input.presence_absence
    models = snakemake.input.models_path
    taxonomy = snakemake.input.taxonomy
    intra_rank_relatedness = snakemake.input.intra_relatedness
    outdir = snakemake.input.outdir

    search_df = read_search(search_result)
    metadata_df = read_metadata(metadata)
    search_df = map_query_to_genome(search_df, metadata_df)
    input_matrix = create_input_matrix(search_df, presence_absence)

    preds = predict_all_ranks(input_matrix, model_base=models, out_dir=outdir)
    print("\n📦 PhagePleats predictions saved to:")
    print(f"    → {os.path.join(outdir, 'taxa_predictions.csv')}")
    preds.to_csv(os.path.join(outdir, 'taxa_predictions.csv'))

    taxonomy_df = pd.read_csv(taxonomy)
    intra_df = pd.read_csv(intra_rank_relatedness)

    dists = compute_closest_training_genomes(presence_absence, input_matrix, taxonomy_df)
    dists.to_csv(os.path.join(outdir, 'phage_closest_hit.csv'), index=False)

    final_summary_df = compute_clade_novelty_summary(
    presence_absence_path=presence_absence,
    input_matrix=input_matrix,
    taxonomy_df=taxonomy_df,
    preds=preds,
    intra_rank_relatedness=intra_df
    )

    final_summary_df.to_csv(os.path.join(outdir, "novel_taxa_summary.csv"))
    print("\n📦 Final novelty-aware prediction summary saved to:")
    print(f"    → {os.path.join(outdir, 'novel_taxa_summary.csv')}")

    def print_happy_phage():
        print(r'''
          ___
         /    \\
        | ^   ^ |     
        |  \_/  |     
         \_____/ 
          |||||       
          |||||       
           |||        
        ___|||___      
       /   |||   \     

    🧬 Predictions are complete. 🎉
        ''')
    print_happy_phage()
