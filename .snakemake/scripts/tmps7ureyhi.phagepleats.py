
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/administrator/programs/miniforge3/envs/phagepleats/lib/python3.10/site-packages', '/home/administrator/.cache/snakemake/snakemake/source-cache/runtime-cache/tmptolzhlcs/file/home/administrator/phd/PhagePleats', '/home/administrator/phd/PhagePleats']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\x06\x07\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8cU/home/administrator/phd/PhagePleats_out_17:04:25_15:07/foldseek_out/search_result.tsv\x94\x8cA/home/administrator/phd/77_test_phages/genome_protein_mapping.csv\x94\x8c\x1cdata/presence_absence.csv.gz\x94\x8c\x0bdata/models\x94\x8c\x11data/taxonomy.csv\x94\x8c\x1fdata/intra_rank_relatedness.csv\x94\x8c6/home/administrator/phd/PhagePleats_out_17:04:25_15:07\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\rsearch_result\x94K\x00N\x86\x94\x8c\x0equery_metadata\x94K\x01N\x86\x94\x8c\x10presence_absence\x94K\x02N\x86\x94\x8c\x0bmodels_path\x94K\x03N\x86\x94\x8c\x08taxonomy\x94K\x04N\x86\x94\x8c\x11intra_relatedness\x94K\x05N\x86\x94\x8c\x06outdir\x94K\x06N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh$\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h*)}\x94\x8c\x05_name\x94h$sNt\x94bh%h(h*\x85\x94R\x94(h*)}\x94h.h%sNt\x94bh\x14h\nh\x16h\x0bh\x18h\x0ch\x1ah\rh\x1ch\x0eh\x1eh\x0fh h\x10ub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8cK/home/administrator/phd/PhagePleats_out_17:04:25_15:07/taxa_predictions.csv\x94\x8cM/home/administrator/phd/PhagePleats_out_17:04:25_15:07/novel_taxa_summary.csv\x94e}\x94(h\x12}\x94(\x8c\ntaxa_preds\x94K\x00N\x86\x94\x8c\x0fclosest_genomes\x94K\x01N\x86\x94uh"]\x94(h$h%eh$h(h*\x85\x94R\x94(h*)}\x94h.h$sNt\x94bh%h(h*\x85\x94R\x94(h*)}\x94h.h%sNt\x94bh<h8h>h9ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x12}\x94h"]\x94(h$h%eh$h(h*\x85\x94R\x94(h*)}\x94h.h$sNt\x94bh%h(h*\x85\x94R\x94(h*)}\x94h.h%sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x12}\x94h"]\x94(h$h%eh$h(h*\x85\x94R\x94(h*)}\x94h.h$sNt\x94bh%h(h*\x85\x94R\x94(h*)}\x94h.h%sNt\x94bub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x04/tmp\x94e}\x94(h\x12}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh"]\x94(h$h%eh$h(h*\x85\x94R\x94(h*)}\x94h.h$sNt\x94bh%h(h*\x85\x94R\x94(h*)}\x94h.h%sNt\x94bhoK\x01hqK\x01hshlub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x12}\x94h"]\x94(h$h%eh$h(h*\x85\x94R\x94(h*)}\x94h.h$sNt\x94bh%h(h*\x85\x94R\x94(h*)}\x94h.h%sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x04PDBs\x94\x8c;/home/administrator/phd/77_test_phages/test_phages_701_pdbs\x94\x8c\x08metadata\x94\x8cA/home/administrator/phd/77_test_phages/genome_protein_mapping.csv\x94\x8c\x06outdir\x94\x8c\x18/home/administrator/phd/\x94u\x8c\x04rule\x94\x8c\x11run_python_script\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c#/home/administrator/phd/PhagePleats\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/home/administrator/phd/PhagePleats/phagepleats.py';
######## snakemake preamble end #########
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

# 4. Create input presence/absence matrix for testing
def create_input_matrix(search, presence_absence_path):
    """
    Creates a presence/absence matrix for input genomes based on detected clusters.

    Args:
        search (pd.DataFrame): Search DataFrame with 'query_genome' and 'target'.
        presence_absence_path (str): Path to presence/absence file with cleaned cluster names.

    Returns:
        pd.DataFrame: Presence/absence matrix (genomes x clusters).
    """
    presence_absence = pd.read_csv(presence_absence_path)
    presence_absence[['accession', 'function']] = presence_absence['cluster_ID_function'].str.split(":", n=1, expand=True)
    presence_absence.set_index('cluster_ID_function', inplace=True)
    cluster_index = presence_absence.index

    genomes = search['query_genome'].dropna().unique()
    input_df = pd.DataFrame(index=cluster_index, columns=genomes)

    for genome in tqdm(genomes, desc="Building presence/absence matrix"):
        genome_clusters = search[search['query_genome'] == genome]['target'].values
        for cluster in genome_clusters:
            match = presence_absence[presence_absence['accession'] == cluster]
            if not match.empty:
                full_id = match.index[0]
                input_df.at[full_id, genome] = 1

    return input_df.astype('float32').fillna(0).T


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
    presence_absence = pd.read_csv(presence_absence_path, index_col=0)
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

# --- Helper: Compute Z-score-based novelty flag ---
def assign_novelty_flags(df, intra_rank_relatedness):
    for rank in ["Order", "Family", "Subfamily", "Genus"]:
        intra = intra_rank_relatedness[intra_rank_relatedness["rank"] == rank][[
            "clade", "intra_avg_shared_proteins", "intra_std_shared_proteins",
            "intra_avg_euclidean", "intra_std_euclidean"
        ]].copy()
        intra["clade"] = intra["clade"].astype(str)
        intra_dict = intra.set_index("clade").to_dict(orient='index')

        z_shared, z_euclid, flag = [], [], []
        for i, row in df.iterrows():
            clade = str(row[rank])
            shared = row.get(f"%_shared_with_predicted_{rank}")
            dist = row.get(f"eucl_dist_to_predicted_{rank}")
            stats = intra_dict.get(clade, {})

            if all(k in stats for k in ["intra_avg_shared_proteins", "intra_std_shared_proteins",
                                        "intra_avg_euclidean", "intra_std_euclidean"]):
                try:
                    z_s = (shared - stats["intra_avg_shared_proteins"]) / stats["intra_std_shared_proteins"]
                    z_d = (dist - stats["intra_avg_euclidean"]) / stats["intra_std_euclidean"]
                except:
                    z_s = z_d = np.nan
            else:
                z_s = z_d = np.nan

            z_shared.append(z_s)
            z_euclid.append(z_d)

            # Assign novelty flag
            if pd.isna(z_s) or pd.isna(z_d):
                flag.append("Unknown")
            elif z_s <= -3 or z_d >= 3:
                flag.append("Potential new order")
            elif z_s <= -2 or z_d >= 2:
                flag.append("Potential new family")
            elif z_s <= -1 or z_d >= 1:
                flag.append("Potential new genus")
            else:
                flag.append("Likely member")

        df[f"{rank}_z_shared"] = z_shared
        df[f"{rank}_z_euclid"] = z_euclid
        df[f"{rank}_novelty_flag"] = flag

    return df

# --- Main novelty-aware function ---
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

    results = Parallel(n_jobs=-1, backend="loky")(
        delayed(process_genome)(genome, input_matrix_f32, input_matrix_bool, training_matrix_bool, clade_members, faiss_indexes, preds)
        for genome in tqdm(input_matrix.index, desc="💫 Calculating per-genome clade similarity")
    )
    end = time.perf_counter()
    print(f"⏱️ Took {end - start:.2f} seconds for intra-clade similarity calculations")

    clade_df = pd.DataFrame(results).set_index("Genome")
    final_df = closest_df.join(clade_df)

    print("🧬 Assigning novelty scores...")
    for rank in ["Order", "Family", "Subfamily", "Genus"]:
        final_df[rank] = preds[rank]

    final_df = assign_novelty_flags(final_df, intra_rank_relatedness)
    final_df.index.name = "Genome"
    final_df = final_df.round(2).sort_index()
    print("✅ Novelty summary complete.")
    return final_df

# #Novelty-aware post-prediction QC layer 👾🧬✨
# def compute_clade_novelty_summary(presence_absence_path, input_matrix, taxonomy_df, preds, intra_rank_relatedness):
#     """
#     Compute novelty-aware summary per input genome based on FAISS distances and Jaccard similarity
#     with predicted taxonomic clades.
#     """

#     print("🔍 Loading training presence/absence matrix...")
#     presence_absence = pd.read_csv(presence_absence_path, index_col=0, compression='gzip')
#     training_matrix = presence_absence.astype('float32').T
#     input_matrix = input_matrix.astype('float32')
#     taxonomy_df = taxonomy_df.set_index("Leaves") 

#     print("🔗 Finding closest training genome (Euclidean distance) to input genome...")
#     index = faiss.IndexFlatL2(training_matrix.shape[1])
#     index.add(training_matrix.values)

#     distances, indices = index.search(input_matrix.values, k=1)
#     closest_df = pd.DataFrame({
#         'Genome': input_matrix.index,
#         'closest_training_genome': training_matrix.index[indices.flatten()],
#         'euclidean_dist_to_closest_hit': np.sqrt(distances.flatten())
#     }).set_index("Genome")

#     print("🔬 Computing % shared proteins + FAISS distance to predicted clades...")
    
#     # Precast matrices
#     training_matrix_f32 = training_matrix.astype('float32')
#     training_matrix_bool = training_matrix.astype(bool)
#     input_matrix_f32 = input_matrix.astype('float32')
#     input_matrix_bool = input_matrix.astype(bool)

#     # Build clade genome index per rank
#     clade_members = defaultdict(dict)
#     for rank in ["Order", "Family", "Subfamily", "Genus"]:
#         for clade, clade_df in taxonomy_df.groupby(rank):
#             matched = clade_df.index.intersection(training_matrix.index)
#             if len(matched) > 0:
#                 clade_members[rank][clade] = matched

#     # Step 3: Build FAISS indexes per clade
#     faiss_indexes = defaultdict(dict)
#     for rank in clade_members:
#         for clade, genomes in clade_members[rank].items():
#             vecs = training_matrix_f32.loc[genomes].values
#             if vecs.shape[0] > 0:
#                 index = faiss.IndexFlatL2(vecs.shape[1])
#                 index.add(vecs)
#                 faiss_indexes[rank][clade] = index

#     def process_genome(genome):
#     row = {"Genome": genome}
#     input_vec = input_matrix_f32.loc[genome].values
#     input_vec_bool = input_matrix_bool.loc[genome].values

#     for rank in ["Order", "Family", "Subfamily", "Genus"]:
#         pred_clade = preds.loc[genome, rank]
#         matched_genomes = clade_members[rank].get(pred_clade, [])

#         if len(matched_genomes) == 0:
#             row[f"%_shared_with_predicted_{rank}"] = np.nan
#             row[f"eucl_dist_to_predicted_{rank}"] = np.nan
#         else:
#             # --- Jaccard (% shared proteins)
#             try:
#                 clade_bools = training_matrix_bool.loc[matched_genomes].values
#                 jaccard_vals = np.mean([
#                     jaccard_score(input_vec_bool, clade_bools[i], average='binary', zero_division=0)
#                     for i in range(clade_bools.shape[0])
#                 ])
#                 row[f"%_shared_with_predicted_{rank}"] = jaccard_vals
#             except:
#                 row[f"%_shared_with_predicted_{rank}"] = np.nan

#             # --- FAISS (Euclidean)
#             try:
#                 faiss_index = faiss_indexes[rank].get(pred_clade)
#                 if faiss_index:
#                     D, _ = faiss_index.search(input_vec.reshape(1, -1), k=faiss_index.ntotal)
#                     row[f"eucl_dist_to_predicted_{rank}"] = np.sqrt(np.mean(D))
#                 else:
#                     row[f"eucl_dist_to_predicted_{rank}"] = np.nan
#             except:
#                 row[f"eucl_dist_to_predicted_{rank}"] = np.nan

#     # results.append(row)
#     start = time.perf_counter()
#     results = []

#     training_matrix_f32 = training_matrix.astype('float32')
#     training_matrix_bool = training_matrix.astype(bool)
#     input_matrix_f32 = input_matrix.astype('float32')
#     input_matrix_bool = input_matrix.astype(bool)

#     clade_members = defaultdict(dict)
#     for rank in clade_members:
#         for clade, genomes in clade_members[rank].items():
#             clade_matrix = training_matrix_f32.loc[genomes].values
#             if clade_matrix.shape[0] > 0:
#                 index = faiss.IndexFlatL2(clade_matrix.shape[1])
#                 index.add(clade_matrix)
#                 faiss_indexes[rank][clade] = index

#     for rank in ["Order", "Family", "Subfamily", "Genus"]:
#         for clade, sub_df in taxonomy_df.groupby(rank):
#             matched = sub_df.index.intersection(training_matrix.index)
#             if len(matched) > 0:
#                 clade_members[rank][clade] = matched

#     for rank in ["Order", "Family", "Subfamily", "Genus"]:
#         preds[rank] = preds[rank].astype(str)
#         taxonomy_df[rank] = taxonomy_df[rank].astype(str)

#     for genome in tqdm(input_matrix.index, desc="💫 Calculating per-genome clade similarity"):
#         row = {"Genome": genome}
#         input_vec = input_matrix.loc[genome].values.astype('float32')

#         for rank in ["Order", "Family", "Subfamily", "Genus"]:
#             pred_clade = preds.loc[genome, rank]
#             genomes_in_clade = taxonomy_df[taxonomy_df[rank] == pred_clade].index
#             matching_genomes = training_matrix.index.intersection(genomes_in_clade)
#             if matching_genomes.empty:
#                 row[f"%_shared_with_predicted_{rank}"] = np.nan
#                 row[f"eucl_dist_to_predicted_{rank}"] = np.nan
#             else:
#                 # % Shared Proteins (Jaccard)
#                 jaccards = [
#                     jaccard_score(
#                         input_vec.astype(bool),
#                         training_matrix.loc[g].values.astype(bool),
#                         average='binary',
#                         zero_division=0
#                     )
#                     for g in matching_genomes
#                 ]
#                 row[f"%_shared_with_predicted_{rank}"] = np.mean(jaccards)

#                 # Euclidean Distance (FAISS)
#                 clade_matrix = training_matrix.loc[matching_genomes].values.astype('float32')
#                 input_vec_2d = input_vec.reshape(1, -1)
#                 clade_index = faiss.IndexFlatL2(clade_matrix.shape[1])
#                 clade_index.add(clade_matrix)
#                 distances, _ = clade_index.search(input_vec_2d, k=clade_matrix.shape[0])
#                 #IndexFlatL2 in FAISS returns L2 squared distances so get the sq root
#                 row[f"eucl_dist_to_predicted_{rank}"] = np.sqrt(np.mean(distances))

#         results.append(row)
#     end = time.perf_counter()
#     print(f"⏱️ Took {end - start:.4} to run intra clade loop....")

#     clade_df = pd.DataFrame(results).set_index("Genome")
#     print("🧬 Merging intra-clade reference statistics...")
#     final_df = closest_df.join(clade_df)
#     original_index = final_df.index.copy()

#     for rank in ["Order", "Family", "Subfamily", "Genus"]:
#         preds[rank] = preds[rank].astype(str)
#         final_df[rank] = preds[rank]

#         intra = intra_rank_relatedness[intra_rank_relatedness["rank"] == rank][
#             ["clade", "intra_avg_shared_proteins", "intra_avg_euclidean"]
#         ].copy()
#         intra["clade"] = intra["clade"].astype(str)

#         intra_dict_shared = dict(zip(intra["clade"], intra["intra_avg_shared_proteins"]))
#         intra_dict_dist = dict(zip(intra["clade"], intra["intra_avg_euclidean"]))

#         final_df[f"{rank}_intra_avg_shared_proteins"] = final_df[rank].map(intra_dict_shared)
#         final_df[f"{rank}_intra_avg_euclidean_dist"] = final_df[rank].map(intra_dict_dist)

#         # Handle NaNs safely in novelty comparisons
#         final_df[f"{rank}_novel_by_shared"] = (
#             final_df[f"%_shared_with_predicted_{rank}"] < final_df[f"{rank}_intra_avg_shared_proteins"]
#         )
#         final_df[f"{rank}_novel_by_distance"] = (
#             final_df[f"eucl_dist_to_predicted_{rank}"] > final_df[f"{rank}_intra_avg_euclidean_dist"]
#         )
 
#     final_df.index = original_index
#     final_df.index.name = "Genome"
#     final_df = final_df.round(2).sort_index()

#     print("✅ Novelty summary complete.")
#     return final_df

# Main Snakemake execution
if __name__ == "__main__":
    def print_phage_ascii():
        print(r'''
          ___
        /     \
       | o   o |
        \  ^  /
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
