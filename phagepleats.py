import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
import faiss

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
        df_proba.to_csv(os.path.join(prob_dir, f"{rank}_proba.csv"), index=False)

    final_df = pd.concat(all_preds, axis=1)
    final_df.index.name = "Genome"
    final_df = final_df.sort_index()

    print(f"\n✅ Final predictions saved to '{os.path.join(out_dir, 'taxonomy_predictions.csv')}'")
    print(f"🧬📂 Full probability matrices saved to '{prob_dir}/'")

    return final_df

# 7. Compute distances to training matrix using FAISS
def compute_closest_training_genomes(presence_absence, input_matrix):
    """
    Finds the closest training genome to each input genome using FAISS (Euclidean distance).

    Args:
        training_matrix_T (pd.DataFrame): Transposed training matrix (genomes x features).
        input_matrix (pd.DataFrame): Transposed input matrix (genomes x features).

    Returns:
        pd.DataFrame: DataFrame containing closest training genome and distance per input genome.
    """
    # Convert to float32 as required by FAISS
    presence_absence = pd.read_csv(presence_absence, index_col=0)
    X_train = presence_absence.astype('float32').T
    X_input = input_matrix.astype('float32')

    # Convert to numpy arrays
    X_train_np = X_train.values
    X_input_np = X_input.values

    # Build FAISS index using Euclidean distance (L2)
    index = faiss.IndexFlatL2(X_train_np.shape[1])
    index.add(X_train_np)

    # Perform nearest neighbor search
    distances, indices = index.search(X_input_np, k=1)

    # Map indices back to genome names
    closest_genomes = X_train.index[indices.flatten()]
    input_ids = X_input.index

    # Construct output DataFrame
    results = pd.DataFrame({
        'input_genome': input_ids,
        'closest_training_genome': closest_genomes,
        'euclidean_distance': distances.flatten()
    })

    return results

#Novelty-aware post-prediction QC layer 👾🧬✨



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
    ''')
        
    # ... call it before prediction starts lol:
    print_phage_ascii()

    search_result = snakemake.input.search_result
    metadata = snakemake.input.query_metadata
    presence_absence = snakemake.input.presence_absence
    models = snakemake.input.models_path
    out_dir = snakemake.input.outdir

    search_df = read_search(search_result)
    metadata_df = read_metadata(metadata)
    search_df = map_query_to_genome(search_df, metadata_df)
    input_matrix = create_input_matrix(search_df, presence_absence)

    preds = predict_all_ranks(input_matrix, model_base=models, out_dir=out_dir)
    preds.to_csv(os.path.join(out_dir, 'taxa_predictions.csv'))

    dists = compute_closest_training_genomes(presence_absence, input_matrix)  # assumes p_a_matrix is defined elsewhere or passed in
    dists.to_csv(os.path.join(out_dir, 'eucl_distances.csv'), index=False)
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
