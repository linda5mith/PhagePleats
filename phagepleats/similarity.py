import pandas as pd
import numpy as np
import faiss
from sklearn.metrics import jaccard_score
from collections import defaultdict


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
    #sqrt_distances = np.sqrt(distances.flatten())

    # Set taxonomy index
    taxonomy_df = taxonomy_df.set_index("Leaves").astype(str)

    # Compute results
    results = []
    for genome, closest in zip(input_matrix.index, closest_genomes):
        print(genome, closest)
        input_vec = input_matrix.loc[genome].astype(bool)
        closest_vec = training_matrix.loc[closest].astype(bool)
        shared = jaccard_score(input_vec, closest_vec, average='binary', zero_division=0)

        taxonomy_df = taxonomy_df.set_index("Leaves").astype(str)
        tax = taxonomy_df.loc[closest].to_dict() if closest in taxonomy_df.index else {}
        print(tax)

        results.append({
            "input_genome": genome,
            "closest_training_genome": closest,
            "%_shared_proteins": shared,
            "Order": tax.get("Order", np.nan),
            "Family": tax.get("Family", np.nan),
            "Subfamily": tax.get("Subfamily", np.nan),
            "Genus": tax.get("Genus", np.nan)
        })

    return pd.DataFrame(results).round(3).sort_values(by='input_genome')

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