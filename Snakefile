configfile: "config.yml"
import os
import yaml
import glob
import json
import pandas as pd
from datetime import datetime
from phagepleats.io import read_search, read_metadata, map_query_to_genome, create_input_matrix
from phagepleats.prediction import predict_all_ranks, hierarchical_lookup
from phagepleats.similarity import compute_closest_training_genomes
from phagepleats.novelty import compute_clade_novelty_summary

# Load config file
with open("config.yml", 'r') as f:
    config = yaml.safe_load(f)

# Get the directory to installation of Snakefile
installation_path = os.path.abspath('.')
input_PDBs = config['PDBs']
cluster_reps = os.path.join(installation_path, 'cluster_reps')
today_date = datetime.now().strftime('%d-%m-%y')
outdir = os.path.join(config['outdir'], f"PhagePleats_out_{today_date}")
log_dir = os.path.join(outdir, "log")
foldseek_out_dir = os.path.join(outdir, "foldseek_out")

# Create necessary directories
os.makedirs(log_dir, exist_ok=True)
os.makedirs(foldseek_out_dir, exist_ok=True)

# Check outdir exists
shell.executable("/bin/bash")
shell.prefix("mkdir -p {outdir} && ")
os.makedirs(log_dir, exist_ok=True)
os.makedirs(foldseek_out_dir, exist_ok=True)

# Print paths for debugging
print(f"Installation Path: {os.path.abspath('.')}")
print(f"Cluster Reps Path: {os.path.join(os.path.abspath('.'), 'cluster_reps')}")
print(f"Output Directory: {outdir}")
print(f"Log Directory: {log_dir}")
print(f"Foldseek Output Directory: {foldseek_out_dir}")

# Defines final targets of workflow
rule all:
    input:
        os.path.join(foldseek_out_dir, 'search_result.tsv'),
        os.path.join(outdir, 'taxa_predictions.csv'),
        os.path.join(outdir, 'taxa_novelty_summary.csv')

# 1. Create foldseek DB of query pdbs
rule foldseek_createDB:
    input:
        path_to_PDBs = input_PDBs
    output:
        queryDB = os.path.join(foldseek_out_dir, 'queryDB')
    log:
        os.path.join(log_dir, "foldseek_createDB.log")
    shell:
        """
        foldseek createdb {input.path_to_PDBs} {output.queryDB} --threads 10 | tee -a {log}
        """

# 2. Foldseek search queryDB against cluster reps targetDB
checkpoint foldseek_search:
    input:
        query=os.path.join(foldseek_out_dir, "queryDB"),
        target="data/targetDB"
    output:
        alignment_dir = directory(os.path.join(outdir, "alignments"))
    log:
        os.path.join(log_dir, "foldseek_search.log")
    shell:
        """
        mkdir -p {output}
        foldseek search {input.query} {input.target} {output}/aln {output}/tmpFolder \
            -c 0.4 --cov-mode 5 --tmscore-threshold 0.3 --max-seqs 5 \
            --exhaustive-search 1 -e 0.01 | tee -a {log}
        """

# 3. Checkpoint to wait for all alignment files to be created
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.foldseek_search.get(**wildcards).output[0]
    alignment_files = glob_wildcards(os.path.join(checkpoint_output, "aln.{i}")).i
    return expand(os.path.join(checkpoint_output, "aln.{i}"), i=alignment_files)

# 4. Convert search alignments to tsv
rule foldseek_convertalis:
    input:
        query=os.path.join(foldseek_out_dir, "queryDB"),
        target="data/targetDB",
        aln_files=aggregate_input
    output:
        result=os.path.join(foldseek_out_dir, "search_result.tsv")
    log:
        os.path.join(log_dir, "foldseek_convertalis.log")
        
    shell:
        """
        foldseek convertalis {input.query} {input.target} {outdir}/alignments/aln {output.result} | tee -a {log}
        """

rule run_phagepleats_pipeline:
    input:
        search_result=os.path.join(foldseek_out_dir, 'search_result.tsv'),
        query_metadata=config['metadata'],
        presence_absence="data/presence_absence.csv.gz",
        models_path="data/models",
        taxonomy_json="data/taxonomy.json",
        taxonomy_df="data/taxonomy.csv",
        intra_relatedness="data/intra_rank_relatedness.csv",
        outdir=outdir
    output:
        taxa_preds=os.path.join(outdir, 'taxa_predictions.csv'),
        closest_genomes=os.path.join(outdir, 'taxa_novelty_summary.csv')
    run:
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
        
        print_phage_ascii()

        # 1: Load and preprocess input data
        search_df = read_search(input.search_result)
        metadata_df = read_metadata(input.query_metadata)
        search_df = map_query_to_genome(search_df, metadata_df)
        input_matrix = create_input_matrix(search_df, input.presence_absence)

        # 2 Load taxonomy
        with open(input.taxonomy_json) as f:
            hierarchy = json.load(f)

        # 3. Predict all ranks
        flat_preds = predict_all_ranks(input_matrix, out_dir=input.outdir)
        preds = hierarchical_lookup(flat_preds, hierarchy)
        preds.to_csv(output.taxa_preds)

        # 4: Compute novelty
        taxonomy_df = pd.read_csv(input.taxonomy_df)
        intra_df = pd.read_csv(input.intra_relatedness)

        novelty_summary_df = compute_clade_novelty_summary(
            input.presence_absence,
            input_matrix=input_matrix,
            taxonomy_df=taxonomy_df,
            preds=preds,
            intra_rank_relatedness=intra_df
        )
        novelty_summary_df.to_csv(output.closest_genomes)

        print("✅ All results saved.")

        def print_happy_phage():
            print(r'''
           ___
         /    \\
        | 0   0 |     
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
        


