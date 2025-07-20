import pandas as pd
import os

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
    if not os.path.exists(path_to_input_metadata):
        raise FileNotFoundError(f"Metadata file not found: {path_to_input_metadata}")

    try:
        metadata = pd.read_csv(path_to_input_metadata)
    except pd.errors.ParserError as e:
        raise ValueError(f"Failed to read metadata as CSV: {e}")

    required_cols = {"protein", "genome"}
    if not required_cols.issubset(metadata.columns):
        raise ValueError(f"Metadata must contain columns: {required_cols}, but found: {set(metadata.columns)}")

    print(f"✅ Loaded input metadata.")
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
    if 'protein' not in metadata.columns or 'genome' not in metadata.columns:
        raise ValueError("Metadata must contain 'protein' and 'genome' columns.")

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
    print("✅ Building presence/absence matrix...")
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

    genomes_with_no_hits = final.index[final.sum(axis=1) == 0].tolist()
    if genomes_with_no_hits:
        print(f"⚠️ {len(genomes_with_no_hits)} input genome(s) had no hits to any cluster:")
        for g in genomes_with_no_hits:
            print(f"   - {g}")
    else:
        print(f"✅ All {len(final)} input genomes have hits to at least one protein cluster.")
    return final    