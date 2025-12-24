import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import pickle
from collections import defaultdict
from importlib.resources import files

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

def build_taxonomy_hierarchy(taxa_df):
    """
    Build a nested taxonomy hierarchy from a DataFrame with Order, Family, Subfamily, and Genus.

    Returns:
        dict: Nested dict representing the hierarchy from Order down to Genus.
    """
    hierarchy = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

    for _, row in taxa_df.iterrows():
        order = row.get('Order')
        family = row.get('Family') if pd.notna(row.get('Family')) else None
        subfamily = row.get('Subfamily') if pd.notna(row.get('Subfamily')) else None
        genus = row.get('Genus') if pd.notna(row.get('Genus')) else None

        if not order or not family or not subfamily or not genus:
            continue

        hierarchy[order]['Family'].add(family)
        if family not in hierarchy[order]:
            hierarchy[order][family] = {'Subfamily': set()}
        hierarchy[order][family]['Subfamily'].add(subfamily)
        if subfamily not in hierarchy[order][family]:
            hierarchy[order][family][subfamily] = {'Genus': set()}
        hierarchy[order][family][subfamily]['Genus'].add(genus)

    return hierarchy

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

    for taxon, model in tqdm(models.items(), desc=f"Predicting {rank}"):
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

def predict_all_ranks(X, ranks=["Order", "Family", "Subfamily", "Genus"], 
                      model_base=None, out_dir="outputs"):
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
    if model_base is None:
        model_base = files("phagepleats") / "resources" / "data" / "models"
        model_base = str(model_base)

    all_preds = []
    os.makedirs(out_dir, exist_ok=True)
    prob_dir = os.path.join(out_dir, "rank_probabilities")
    os.makedirs(prob_dir, exist_ok=True)

    for rank in ranks:
        # print(f"\nPredicting for rank: {rank}")
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


def hierarchical_lookup(preds_df, hierarchy):
    """
    For each genome, traverse the predicted hierarchy and retain only valid paths.
    """
    cleaned_preds = []

    for genome, row in preds_df.iterrows():
        result = {"Genome": genome}

        order = row.get("Order")
        if order not in hierarchy:
            result.update({rank: None for rank in ["Order", "Order_prob", "Family", "Family_prob", "Subfamily", "Subfamily_prob", "Genus", "Genus_prob"]})
            cleaned_preds.append(result)
            continue

        result["Order"] = order
        result["Order_prob"] = row.get("Order_prob")

        # ----------------------
        # FAMILY
        # ----------------------
        possible_families = hierarchy[order].keys()
        family = row.get("Family") if row.get("Family") in possible_families else None
        result["Family"] = family
        result["Family_prob"] = row.get("Family_prob") if family else None

        # ----------------------
        # SUBFAMILY
        # ----------------------
        subfamily = None
        if family:
            possible_subfamilies = hierarchy[order][family].keys()
            subfamily = row.get("Subfamily") if row.get("Subfamily") in possible_subfamilies else None
        result["Subfamily"] = subfamily
        result["Subfamily_prob"] = row.get("Subfamily_prob") if subfamily else None

        # ----------------------
        # GENUS
        # ----------------------
        genus = None
        if subfamily:
            possible_genera = hierarchy[order][family][subfamily]
            genus = row.get("Genus") if row.get("Genus") in possible_genera else None
        result["Genus"] = genus
        result["Genus_prob"] = row.get("Genus_prob") if genus else None

        cleaned_preds.append(result)

    return pd.DataFrame(cleaned_preds).set_index("Genome")

