import pandas as pd
import numpy as np
import time
from tqdm import tqdm
from collections import defaultdict
import faiss
from .similarity import preprocess_matrices, build_faiss_indexes, process_genome

def classify_novelty(z_s, z_d, rank):
    if pd.isna(z_s) or pd.isna(z_d):
        return "Unknown"
    for threshold, label in zip([3, 2, 1], ["order", "family", "subfamily", "genus"]):
        if z_s <= -threshold or z_d >= threshold:
            return f"Potential new {label if rank.lower() == label else rank.lower()}"
    return "Likely member"

def compute_z_scores_and_flag_all_ranks(phage_df, intra_df):
    for rank in ["Order", "Family", "Subfamily", "Genus"]:
        z_shared_list = []
        flag_list = []

        for _, row in phage_df.iterrows():
            predicted_clade = row[f"Closest_{rank}"]
            clade_metrics = intra_df[(intra_df['rank'] == rank) & (intra_df['clade'] == predicted_clade)]

            if clade_metrics.empty:
                z_shared = None
                novelty_flag = "Clade not in intra table"
            else:
                mean_shared = clade_metrics["intra_avg_shared_proteins"].values[0]
                std_shared = clade_metrics["intra_std_shared_proteins"].values[0]

                try:
                    shared = row[f"%_shared_with_predicted_{rank}"]
                    z_shared = (shared - mean_shared) / std_shared if std_shared != 0 else np.nan
                except:
                    z_shared = np.nan

                novelty_flag = classify_novelty(z_shared, np.nan, rank)

            z_shared_list.append(z_shared)
            flag_list.append(novelty_flag)

        phage_df[f"{rank}_z_shared_proteins"] = z_shared_list
        phage_df[f"{rank}_novelty_flag"] = flag_list

    return phage_df


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

    # Add closest genome's taxonomy to df
    closest_tax = taxonomy_df.loc[closest_df['closest_training_genome']]
    for rank in ["Order", "Family", "Subfamily", "Genus"]:
        closest_df[f"Closest_{rank}"] = closest_tax[rank].values

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

    # 🧹 Keep only essential columns
    cols_to_keep = [
        "closest_training_genome",
        "Closest_Order", "Closest_Family", "Closest_Subfamily", "Closest_Genus",
        "%_shared_with_predicted_Order", "%_shared_with_predicted_Family",
        "%_shared_with_predicted_Subfamily", "%_shared_with_predicted_Genus",
        "Order_z_shared_proteins", "Order_novelty_flag",
        "Family_z_shared_proteins", "Family_novelty_flag",
        "Subfamily_z_shared_proteins", "Subfamily_novelty_flag",
        "Genus_z_shared_proteins", "Genus_novelty_flag"
    ]
    final_df = final_df[cols_to_keep]

    print("✅ Novelty summary complete.")
    return final_df
