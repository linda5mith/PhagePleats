import pandas as pd
import numpy as np
import time
from tqdm import tqdm
from collections import defaultdict
from .similarity import preprocess_matrices, build_faiss_indexes, process_genome
import faiss

def classify_novelty(z_s, rank):
    """
    Classify the novelty of a genome at a given taxonomic rank based on z-scores.

    A genome is considered potentially novel if its z-score for % shared proteins
    is substantially below the intra-clade average (i.e., more than 1, 2, or 3 standard deviations).

    Parameters
    ----------
    z_s : float or np.nan
        The z-score of the % shared proteins between the input genome and the predicted clade.
        Lower values indicate less similarity and higher novelty.
    rank : str
        The taxonomic rank being evaluated (e.g., "Order", "Family", "Subfamily", "Genus").

    Returns
    -------
    str
        A string label indicating the novelty classification:
        - "Potential new order/family/subfamily/genus – highly/very/slightly novel" for z_s ≤ -3, -2, or -1
        - "Likely member" if z_s is within expected range i.e. -1 <= z_s <= 1
        - "Unknown" if z_s is NaN
    """
    if pd.isna(z_s):
        return "Unknown"
    if z_s <= -3:
        return f"Potential new {rank.lower()} – highly novel"
    elif z_s <= -2:
        return f"Potential new {rank.lower()} – very novel"
    elif z_s <= -1:
        return f"Potential new {rank.lower()} – slightly novel"
    else:
        return f"Likely {rank.lower()} member"


def compute_z_scores_and_flag_all_ranks(phage_df, intra_df):
    """
    Compute z-scores for % shared proteins and assign novelty flags for each taxonomic rank.

    For each genome and taxonomic rank (Order, Family, Subfamily, Genus), this function:
    - Retrieves the predicted clade.
    - Looks up the intra-clade mean and standard deviation of % shared proteins.
    - Calculates the z-score for the genome's % shared proteins with its predicted clade.
    - Uses the z-score to assign a novelty flag (e.g., 'Likely member', 'Potential new genus').
    - Adds two new columns per rank: one for the z-score and one for the novelty classification.

    Parameters
    ----------
    phage_df : pandas.DataFrame
        A DataFrame containing genome-level predictions, including:
        - 'Closest_<Rank>' columns for each taxonomic rank
        - '%_shared_with_predicted_<Rank>' columns for % shared protein similarity

    intra_df : pandas.DataFrame
        A reference table of intra-clade metrics containing:
        - 'rank': taxonomic rank name
        - 'clade': clade name
        - 'intra_avg_shared_proteins': mean % shared proteins within the clade
        - 'intra_std_shared_proteins': standard deviation of % shared proteins within the clade

    Returns
    -------
    phage_df : pandas.DataFrame
        The input DataFrame with eight additional columns:
        - '<Rank>_z_shared_proteins': z-score of % shared proteins for the predicted clade
        - '<Rank>_novelty_flag': novelty classification label based on the z-score
    """
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

                novelty_flag = classify_novelty(z_shared, rank)

            z_shared_list.append(z_shared)
            flag_list.append(novelty_flag)

        phage_df[f"{rank}_shared_proteins_Z_score"] = z_shared_list
        phage_df[f"{rank}_novelty_flag"] = flag_list

    return phage_df


def compute_clade_novelty_summary(presence_absence_path, input_matrix, taxonomy_df, preds, intra_rank_relatedness):
    """
    Generate a novelty summary for input genomes based on shared protein content and clade-level comparisons.

    This function evaluates how similar each input genome is to known clades at multiple taxonomic ranks 
    (Order, Family, Subfamily, Genus) by:
    - Finding the closest training genome by Euclidean distance in a presence/absence matrix.
    - Assigning predicted taxonomy from trained classifiers.
    - Calculating % shared proteins and distance to predicted clades using FAISS.
    - Computing z-scores for each genome’s similarity to its predicted clade.
    - Assigning novelty flags based on these z-scores.

    Parameters
    ----------
    presence_absence_path : str
        Path to the gzipped CSV file containing the training presence/absence matrix
        (rows = protein clusters, columns = genomes).

    input_matrix : pandas.DataFrame
        Presence/absence matrix of input genomes to evaluate (rows = genomes, columns = protein clusters).

    taxonomy_df : pandas.DataFrame
        DataFrame with taxonomic annotations for training genomes.
        Must include a 'Leaves' column (genome IDs) and taxonomic ranks: 'Order', 'Family', 'Subfamily', 'Genus'.

    preds : dict
        Dictionary mapping taxonomic ranks to predicted clades for each input genome.
        Keys: 'Order', 'Family', 'Subfamily', 'Genus'.
        Values: pandas Series or dict-like objects indexed by genome name.

    intra_rank_relatedness : pandas.DataFrame
        DataFrame containing intra-clade similarity statistics:
        - Columns: 'rank', 'clade', 'intra_avg_shared_proteins', 'intra_std_shared_proteins'

    Returns
    -------
    final_df : pandas.DataFrame
        A summary DataFrame indexed by genome, with the following columns:
        - closest_training_genome: Nearest neighbor genome by Euclidean distance
        - Closest_<Rank>: Closest clade by taxonomy from nearest neighbor
        - %_shared_with_predicted_<Rank>: % of proteins shared with predicted clade
        - <Rank>_shared_proteins_Z_score: Z-score of shared protein similarity
        - <Rank>_novelty_flag: Novelty classification label (e.g., 'Potential new genus – highly novel')
    """
    print("🔍 Loading training presence/absence matrix...")
    presence_absence = pd.read_csv(presence_absence_path, index_col=0, compression='gzip')
    training_matrix = presence_absence.astype('float32').T
    input_matrix = input_matrix.astype('float32')
    taxonomy_df = taxonomy_df.set_index("Leaves")

    print("🔗 Finding closest training genome to input genome...")
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
    training_matrix_f32, training_matrix_bool, input_matrix_f32, input_matrix_bool, clade_members = preprocess_matrices(
        training_matrix, input_matrix, taxonomy_df, preds
    )
    faiss_indexes = build_faiss_indexes(clade_members, training_matrix_f32)

    results = []
    for genome in tqdm(input_matrix.index, desc="💫 Calculating per-genome clade similarity"):
        row = process_genome(genome, input_matrix_f32, input_matrix_bool, training_matrix_bool, clade_members, faiss_indexes, preds)
        results.append(row)

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
        "Order_shared_proteins_Z_score", "Order_novelty_flag",
        "Family_shared_proteins_Z_score", "Family_novelty_flag",
        "Subfamily_shared_proteins_Z_score", "Subfamily_novelty_flag",
        "Genus_shared_proteins_Z_score", "Genus_novelty_flag"
    ]
    final_df = final_df[cols_to_keep]

    print("✅ Novelty summary complete.")
    return final_df
