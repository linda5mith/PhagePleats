import pandas as pd
import os

############################################
# 1. Read Foldseek result file
############################################

def read_search(path_to_foldseek_search):
    """
    Reads Foldseek search results and removes '.pdb' extensions
    from query and target columns.
    """
    search = pd.read_csv(path_to_foldseek_search, sep="\t", header=None)
    search.columns = [
        "query", "target", "fident", "alnlen", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend", "evalue", "bits"
    ]
    search["query"] = search["query"].str.replace(".pdb", "", regex=True)
    search["target"] = search["target"].str.replace(".pdb", "", regex=True)
    return search


############################################
# 2. Read query protein → genome metadata
############################################

def read_metadata(path_to_input_metadata):
    """
    Reads metadata mapping query proteins to query genomes.
    """
    if not os.path.exists(path_to_input_metadata):
        raise FileNotFoundError(f"Metadata file not found: {path_to_input_metadata}")

    metadata = pd.read_csv(path_to_input_metadata)

    required = {"protein", "genome"}
    if not required.issubset(metadata.columns):
        raise ValueError(
            f"Metadata must contain columns {required}, "
            f"but found {set(metadata.columns)}"
        )

    print("Loaded input metadata.")
    return metadata


############################################
# 3. Read genome-level taxonomy
############################################

def read_taxonomy():
    """
    Loads built-in genome-level taxonomy table.
    """
    taxonomy_path = os.path.join(
        os.path.dirname(__file__),
        "resources",
        "data",
        "taxonomy.csv",
    )

    if not os.path.exists(taxonomy_path):
        raise FileNotFoundError(f"Taxonomy file not found: {taxonomy_path}")

    taxonomy = pd.read_csv(taxonomy_path)

    if "Leaves" not in taxonomy.columns:
        raise ValueError(
            "Taxonomy file must contain a 'Leaves' column (genome accession)."
        )

    return taxonomy


############################################
# 4. Map query proteins → query genomes
############################################

def map_query_to_genome(search, metadata):
    """
    Adds a 'query_genome' column to the Foldseek hits.
    """
    prot2genome = dict(zip(metadata["protein"], metadata["genome"]))
    search = search.copy()
    search["query_genome"] = search["query"].map(prot2genome)
    return search


############################################
# 5. Add cluster rep metadata + target genome
############################################

def add_cluster_protein_metadata(hits, cluster_meta):
    """
    Adds cluster representative metadata AND reference genome
    to Foldseek hit-level table.
    """
    required = {"protein", "cluster_genome_accn"}
    if not required.issubset(cluster_meta.columns):
        raise ValueError(
            f"cluster metadata must contain columns {required}"
        )

    enriched = hits.merge(
        cluster_meta,
        how="left",
        left_on="target",
        right_on="protein",
    )

    # Reference genome for aggregation
    enriched["target_genome"] = enriched["cluster_genome_accn"]
    enriched = enriched.sort_values(by='query_genome')

    return enriched


############################################
# 6. Aggregate hits → genome-level similarity
############################################

def aggregate_hits_by_genome(hits, query_metadata):
    """
    Aggregates protein-level hits into genome-level similarity metrics
    and computes TRUE query-normalized structural coverage.
    """
    hits = hits.copy()

    # Per-alignment query coverage (QC metric)
    hits["query_aln_cov"] = (
        (hits["qend"] - hits["qstart"] + 1) / hits["alnlen"]
    )

    grouped = (
        hits
        .groupby(["query_genome", "target_genome"])
        .agg(
            n_hits=("target", "count"),
            n_query_proteins=("query", "nunique"),
            mean_identity=("fident", "mean"),
            median_identity=("fident", "median"),
            mean_bits=("bits", "mean"),
            sum_bits=("bits", "sum"),
            mean_aln_len=("alnlen", "mean"),
            total_aln_len=("alnlen", "sum"),
            mean_query_cov=("query_aln_cov", "mean"),
        )
        .reset_index()
    )

    # Total query proteins per genome
    query_totals = (
        query_metadata
        .groupby("genome")["protein"]
        .nunique()
        .rename("total_query_proteins")
    )

    grouped = grouped.merge(
        query_totals,
        left_on="query_genome",
        right_index=True,
        how="left",
    )

    # Structural coverage (query-normalized)
    grouped["structural_coverage_pct"] = (
        100 * grouped["n_query_proteins"] / grouped["total_query_proteins"]
    ).round(1)

    return grouped


############################################
# 7. Add taxonomy to genome summary
############################################

def add_taxonomy_to_genome_summary(genome_summary, taxonomy):
    """
    Joins genome-level taxonomy onto aggregated similarity table.
    """
    enriched = genome_summary.merge(
        taxonomy[["Leaves", "Order", "Family", "Subfamily", "Genus"]],
        how="left",
        left_on="target_genome",
        right_on="Leaves",
    )
    return enriched


############################################
# 8. Create presence/absence matrix
############################################

def create_input_matrix(search, presence_absence_path):
    """
    Creates a presence/absence matrix for input genomes
    based on detected clusters.
    """
    print("Building presence/absence matrix...")

    presence_absence = pd.read_csv(presence_absence_path)
    presence_absence[["accession", "function"]] = (
        presence_absence["cluster_ID_function"]
        .str.split(":", n=1, expand=True)
    )
    presence_absence.set_index("cluster_ID_function", inplace=True)

    accession_to_cluster = presence_absence["accession"].to_dict()
    cluster_lookup = {v: k for k, v in accession_to_cluster.items()}

    search = search.copy()
    search["matched_cluster"] = search["target"].map(cluster_lookup)

    valid = search.dropna(subset=["query_genome", "matched_cluster"])

    pivot = (
        valid
        .drop_duplicates(subset=["query_genome", "matched_cluster"])
        .assign(present=1)
        .pivot(
            index="query_genome",
            columns="matched_cluster",
            values="present",
        )
        .fillna(0)
        .astype("float32")
    )

    final = pivot.reindex(
        columns=presence_absence.index,
        fill_value=0,
    )

    genomes_with_no_hits = final.index[final.sum(axis=1) == 0].tolist()
    if genomes_with_no_hits:
        print(
            f"⚠️ {len(genomes_with_no_hits)} input genome(s) had no hits:"
        )
        for g in genomes_with_no_hits:
            print(f"   - {g}")
    else:
        print(
            f"✅ All {len(final)} input genomes have hits to at least one cluster."
        )

    return final
