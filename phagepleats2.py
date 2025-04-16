import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
import faiss
from sklearn.metrics import jaccard_score

RANKS = ["Order", "Family", "Subfamily", "Genus"]

# 1. Read Foldseek result file
def read_search(path):
    search = pd.read_csv(path, sep='\t', header=None, names=[
        'query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen',
        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits'
    ])
    search['query'] = search['query'].str.replace('.pdb', '', regex=True)
    search['target'] = search['target'].str.replace('.pdb', '', regex=True)
    return search

# 2. Read metadata mapping protein -> genome
def read_metadata(path):
    return pd.read_csv(path)

# 3. Map each query to its genome
def map_query_to_genome(search, metadata):
    mapping = dict(zip(metadata['protein'], metadata['genome']))
    search['query_genome'] = search['query'].map(mapping)
    return search

# 4. Create input presence/absence matrix
def create_input_matrix(search, pa_path):
    pa = pd.read_csv(pa_path)
    pa[['accession', 'function']] = pa['cluster_ID_function'].str.split(":", n=1, expand=True)
    pa.set_index('cluster_ID_function', inplace=True)

    genomes = search['query_genome'].dropna().unique()
    input_df = pd.DataFrame(0, index=pa.index, columns=genomes, dtype='float32')

    for genome in tqdm(genomes, desc="Building presence/absence matrix"):
        clusters = search.loc[search['query_genome'] == genome, 'target']
        matched = pa[pa['accession'].isin(clusters)].index
        input_df.loc[matched, genome] = 1

    return input_df.T

# 5. Load and use prediction models
def load_models(path):
    return {
        fname.replace(".pkl", ""): pickle.load(open(os.path.join(path, fname), "rb"))
        for fname in os.listdir(path) if fname.endswith(".pkl")
    }

def predict_from_models(X, models, rank):
    probs = [model.predict_proba(X)[:, 1] for model in models.values()]
    prob_matrix = pd.DataFrame(np.column_stack(probs), index=X.index, columns=models.keys())
    pred_labels = prob_matrix.idxmax(axis=1)
    pred_scores = prob_matrix.max(axis=1)

    df_out = pd.DataFrame({
        f"{rank}": pred_labels,
        f"{rank}_prob": pred_scores
    }, index=X.index)

    return df_out.sort_index(), prob_matrix.sort_index()

def predict_all_ranks(X, ranks=RANKS, model_base="data/models", out_dir="outputs"):
    os.makedirs(out_dir, exist_ok=True)
    prob_dir = os.path.join(out_dir, "rank_probabilities")
    os.makedirs(prob_dir, exist_ok=True)

    all_preds = []
    for rank in ranks:
        print(f"\n🔮🧬✨ Predicting for rank: {rank}")
        path = os.path.join(model_base, rank)
        if not os.path.exists(path):
            print(f"❌ No model directory for {rank}: {path}")
            continue

        models = load_models(path)
        if not models:
            print(f"⚠️ No models found in {path}")
            continue

        df_rank, df_proba = predict_from_models(X, models, rank)
        all_preds.append(df_rank)
        df_proba.to_csv(os.path.join(prob_dir, f"{rank}_proba.csv"))

    final_df = pd.concat(all_preds, axis=1)
    final_df.index.name = "Genome"
    return final_df.sort_index()

# 6. Closest genome using FAISS
def compute_closest_training_genomes(pa_path, input_matrix):
    pa = pd.read_csv(pa_path, index_col=0).astype('float32').T
    X_train, X_input = pa.values, input_matrix.astype('float32').values

    index = faiss.IndexFlatL2(X_train.shape[1])
    index.add(X_train)
    dists, idxs = index.search(X_input, k=1)

    return pd.DataFrame({
        'input_genome': input_matrix.index,
        'closest_training_genome': pa.index[idxs.flatten()],
        'euclidean_distance': np.sqrt(dists.flatten())
    })

# 7. Novelty-aware post-prediction QC layer
def compute_clade_novelty_summary(pa_path, input_matrix, taxonomy_df, preds, intra_df):
    print("🔍 Loading training presence/absence matrix...")
    pa = pd.read_csv(pa_path, index_col=0).astype('float32').T
    input_matrix = input_matrix.astype('float32')
    taxonomy_df = taxonomy_df.set_index("Leaves")

    print("🔗 Finding closest training genome...")
    index = faiss.IndexFlatL2(pa.shape[1])
    index.add(pa.values)
    dists, idxs = index.search(input_matrix.values, k=1)
    closest_df = pd.DataFrame({
        'Genome': input_matrix.index,
        'closest_training_genome': pa.index[idxs.flatten()],
        'euclidean_dist_to_closest_hit': np.sqrt(dists.flatten())
    }).set_index("Genome")

    print("🔬 Computing % shared proteins + FAISS distance to predicted clades...")
    results = []
    for genome in tqdm(input_matrix.index, desc="💫 Calculating per-genome clade similarity"):
        row = {"Genome": genome}
        input_vec = input_matrix.loc[genome].values.astype('float32')

        for rank in RANKS:
            pred = preds.loc[genome, rank]
            clade_genomes = taxonomy_df[taxonomy_df[rank] == pred].index
            matched = pa.index.intersection(clade_genomes)

            if matched.empty:
                row[f"%_shared_with_predicted_{rank}"] = np.nan
                row[f"eucl_dist_to_predicted_{rank}"] = np.nan
                continue

            jaccards = [jaccard_score(input_vec.astype(bool), pa.loc[g].values.astype(bool), average='binary', zero_division=0) for g in matched]
            row[f"%_shared_with_predicted_{rank}"] = np.mean(jaccards)

            clade_matrix = pa.loc[matched].values.astype('float32')
            faiss_index = faiss.IndexFlatL2(clade_matrix.shape[1])
            faiss_index.add(clade_matrix)
            dist, _ = faiss_index.search(input_vec.reshape(1, -1), k=clade_matrix.shape[0])
            row[f"eucl_dist_to_predicted_{rank}"] = np.sqrt(np.mean(dist))

        results.append(row)

    clade_df = pd.DataFrame(results).set_index("Genome")
    final_df = closest_df.join(clade_df)

    print("🧬 Merging intra-clade reference statistics...")
    for rank in RANKS:
        final_df[rank] = preds[rank].astype(str)
        intra = intra_df[intra_df["rank"] == rank]
        intra_dict_shared = dict(zip(intra["clade"], intra["intra_avg_shared_proteins"]))
        intra_dict_dist = dict(zip(intra["clade"], intra["intra_avg_euclidean"]))

        final_df[f"{rank}_intra_avg_shared_proteins"] = final_df[rank].map(intra_dict_shared)
        final_df[f"{rank}_intra_avg_euclidean_dist"] = final_df[rank].map(intra_dict_dist)

        final_df[f"{rank}_novel_by_shared"] = final_df[f"%_shared_with_predicted_{rank}"] < final_df[f"{rank}_intra_avg_shared_proteins"]
        final_df[f"{rank}_novel_by_distance"] = final_df[f"eucl_dist_to_predicted_{rank}"] > final_df[f"{rank}_intra_avg_euclidean_dist"]

    return final_df.round(2).sort_index()

# CLI Entrypoint
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

    print_phage_ascii()

    search_result = snakemake.input.search_result
    metadata = snakemake.input.query_metadata
    pa_path = snakemake.input.presence_absence
    models_path = snakemake.input.models_path
    taxonomy_path = snakemake.input.taxonomy
    intra_path = snakemake.input.intra_relatedness
    outdir = snakemake.input.outdir

    search_df = read_search(search_result)
    metadata_df = read_metadata(metadata)
    search_df = map_query_to_genome(search_df, metadata_df)
    input_matrix = create_input_matrix(search_df, pa_path)

    preds = predict_all_ranks(input_matrix, model_base=models_path, out_dir=outdir)
    preds.to_csv(os.path.join(outdir, 'taxa_predictions.csv'))

    taxonomy_df = pd.read_csv(taxonomy_path)
    intra_df = pd.read_csv(intra_path)

    final_df = compute_clade_novelty_summary(pa_path, input_matrix, taxonomy_df, preds, intra_df)
    final_df.to_csv(os.path.join(outdir, "novel_taxa_summary.csv"))

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

    \U0001f9ec Predictions are complete. 🎉
        ''')

    print("\n📦 Final novelty-aware prediction summary saved to:")
    print(f"    → {os.path.join(outdir, 'novel_taxa_summary.csv')}")
    print_happy_phage()
