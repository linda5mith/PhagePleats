
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/administrator/programs/miniforge3/envs/CaudoAuto/lib/python3.10/site-packages', '/home/administrator/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpwumabvpj/file/home/administrator/phd/PhagePleats', '/home/administrator/phd/PhagePleats']); import pickle; snakemake = pickle.loads(b"\x80\x04\x95k\x05\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8cF/home/administrator/phd/PhagePleats_out/foldseek_out/search_result.tsv\x94\x8c)data/presence_absence_c_40_tm_30_4100.csv\x94\x8c\x0bdata/models\x94\x8c'/home/administrator/phd/PhagePleats_out\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\rsearch_result\x94K\x00N\x86\x94\x8c\x10presence_absence\x94K\x01N\x86\x94\x8c\x0bmodels_path\x94K\x02N\x86\x94\x8c\x06outdir\x94K\x03N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x1b\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h!)}\x94\x8c\x05_name\x94h\x1bsNt\x94bh\x1ch\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1csNt\x94bh\x11h\nh\x13h\x0bh\x15h\x0ch\x17h\rub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8c</home/administrator/phd/PhagePleats_out/taxa_predictions.csv\x94a}\x94(h\x0f}\x94\x8c\ntaxa_preds\x94K\x00N\x86\x94sh\x19]\x94(h\x1bh\x1ceh\x1bh\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1bsNt\x94bh\x1ch\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1csNt\x94bh2h/ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x0f}\x94h\x19]\x94(h\x1bh\x1ceh\x1bh\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1bsNt\x94bh\x1ch\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1csNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0f}\x94h\x19]\x94(h\x1bh\x1ceh\x1bh\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1bsNt\x94bh\x1ch\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1csNt\x94bub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x04/tmp\x94e}\x94(h\x0f}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x19]\x94(h\x1bh\x1ceh\x1bh\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1bsNt\x94bh\x1ch\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1csNt\x94bhcK\x01heK\x01hgh`ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0f}\x94h\x19]\x94(h\x1bh\x1ceh\x1bh\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1bsNt\x94bh\x1ch\x1fh!\x85\x94R\x94(h!)}\x94h%h\x1csNt\x94bub\x8c\x06config\x94}\x94(\x8c\x04PDBs\x94\x8cA/home/administrator/phd/soil_2024/circular_irish_L001_2_3_4_folds\x94\x8c\x06outdir\x94\x8c\x18/home/administrator/phd/\x94u\x8c\x04rule\x94\x8c\x11run_python_script\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c#/home/administrator/phd/PhagePleats\x94ub."); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/home/administrator/phd/PhagePleats/phagepleats.py';
######## snakemake preamble end #########
from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import re
import argparse
import subprocess
from tqdm import tqdm
from io import StringIO
from rich.console import Console
from rich.progress import track, Progress
from rich.columns import Columns
from rich.markdown import Markdown
from rich.panel import Panel
from datetime import datetime
from sklearn.multioutput import MultiOutputClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
import pickle
import gzip
pd.set_option("future.no_silent_downcasting", True)

#1 Read in result_cluster.tsv
def read_search(path_to_foldseek_search):
    search = pd.read_csv(path_to_foldseek_search,
                 sep='\t',header=None)
    search.columns = ['query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits']
    search['query'] = search['query'].str.replace('.pdb', '', regex=True)
    search['target'] = search['target'].str.replace('.pdb', '', regex=True)
    search.loc[:, 'query'] = search['query'].str.replace('.pdb', '', regex=True)
    search.loc[:,'target'] = search['target'].str.replace('.pdb', '', regex=True)
    search.loc[:, 'query_genome'] = search['query'].str.rsplit('_', n=1).str[0]
    print(search.head())
    return search

#2 Create presence/absence table from inputted genomes
def create_input_matrix(search, p_a_matrix):
    presence_absence = pd.read_csv(p_a_matrix)
    index = presence_absence['cluster_ID_function']
    unique_genomes = search['query_genome'].unique()
    input_df = pd.DataFrame(index=index, columns=unique_genomes)
    cluster_reps = search['target'].unique()
    genomes = search['query_genome'].unique()
    # Loop through clusters and update presence DataFrame
    for genome in tqdm(genomes, desc="Checking genome presence"):
        genome_clusters = search[search['query_genome'] == genome]['target']
        if not genome_clusters.empty:
            # Create a mask for the index
            mask = input_df.index.str.contains('|'.join(genome_clusters), na=False)
            # Update the presence DataFrame
            input_df.loc[mask, genome] = 1
    input_df_T = input_df.fillna(0).T
    return input_df_T

def predict_taxa(input_df, path_to_models):
    output_df = input_df.reset_index()[['index']].rename(columns={'index': 'Genome'})
    model_files = [f for f in os.listdir(path_to_models) if f.endswith('.gz')]
    print(model_files)
    model_path = os.path.join(path_to_models, model_files[0])
    print(model_path)
    with gzip.open(model_path, 'rb') as f:
        print('Loading model...')
        model = pickle.load(f)
    print(f'Predicting with model: {input_df}')
    new_taxa_predictions = model.predict(input_df)
    taxa_df = pd.DataFrame(new_taxa_predictions, columns=['Order', 'Family', 'Genus', 'Species'])
    preds_df = pd.concat([output_df, taxa_df], axis=1)
    return preds_df

#4 Calculate euclidean distances between input genomes and base 4100 matrix
def get_distances(input_df, p_a_matrix):
    print('reaching here....')
    features = pd.read_csv(p_a_matrix, index_col='cluster_ID_function')
    features = features.T

    # Perform PCA
    pca = PCA(n_components=2000)
    reduced_matrix = pca.fit_transform(features)
    
    # Convert reduced_matrix to DataFrame
    reduced_matrix_df = pd.DataFrame(reduced_matrix, index=features.index, columns=[f'PC{i+1}' for i in range(2000)])
    
    # Transform new genome features
    new_genome_features = input_df
    new_genome_reduced = pca.transform(new_genome_features)

    # Convert new_genome_reduced to DataFrame
    new_reduced_df = pd.DataFrame(new_genome_reduced, index=new_genome_features.index, columns=[f'PC{i+1}' for i in range(2000)])
    
    # Print DataFrame heads and shapes
    print(new_reduced_df.head())
    print(new_reduced_df.shape)
    print(reduced_matrix_df.head())
    print(reduced_matrix_df.shape)
    
    # Concatenate DataFrames
    combined_reduced_df = pd.concat([new_reduced_df, reduced_matrix_df])
    
    # Calculate Euclidean distance matrix
    distance_matrix = cdist(combined_reduced_df, combined_reduced_df, metric='euclidean')
    distance_matrix_df = pd.DataFrame(distance_matrix, index=combined_reduced_df.index, columns=combined_reduced_df.index)
    distance_matrix_df = distance_matrix_df.round(4)
    return distance_matrix_df

#5. Convert to PHYLIP
def phylip(dist_matrix):
    matrix = dist_matrix.pivot_table(
            columns='qseqid_genome_ID', 
            index='sseqid_genome_ID', 
            values='distance'
        ).reset_index()
    matrix = matrix.set_index('sseqid_genome_ID')
    pbar.update(1)
    
    # Check for and handle missing rows and columns
    cols = set(matrix.columns.to_list())
    indx = set(matrix.index.to_list())
    cols_missing = list(set(indx) - set(cols))
    idx_missing = list(set(cols) - set(indx))
    
    # Add missing rows
    new_row = pd.Series(np.full(shape=len(idx_missing), fill_value=np.nan), index=idx_missing)
    indexes_insert = pd.DataFrame(index=new_row.index)
    matrix = pd.concat([matrix, indexes_insert], axis=1)
    pbar.update(1)
    # Add missing columns
    missing_df = pd.DataFrame(1, columns=cols_missing, index=matrix.index)
    matrix = pd.concat([matrix, missing_df], axis=1)
    sorted_index = sorted(matrix.index)
    matrix = matrix.loc[sorted_index]
    matrix = matrix[sorted_index]
    matrix.columns = [''] * len(matrix.columns)
    matrix.index.names = [len(matrix)]
    matrix.index = matrix.index.str.replace(' ', '_')
    matrix.index = matrix.index.str.replace(r'\,|\(|\)|\:', '_', regex=True)
    pbar.update(1)
    
    # Ensure taxon names are not longer than 64 characters
    matrix.index = matrix.index.str[:64]
    np.fill_diagonal(matrix.values, 0)
    matrix = matrix.fillna(1).astype(float).round(5).clip(lower=0)
    
    # Save the matrix to a PHYLIP file
    phylip_file = os.path.join(self.out_dir, 'phylip.dist')
    matrix.to_csv(phylip_file, header=True, index=True, sep=' ')
    return matrix
    
#6. Decenttree

#7: % of shared proteins between new genome and existing ones

if __name__ == "__main__":
    # Access the input file
    search_result = snakemake.input.search_result
    search_result = read_search(search_result)
    p_a_matrix = snakemake.input.presence_absence
    input_df = create_input_matrix(search_result, p_a_matrix)
    models = snakemake.input.models_path
    out_dir = snakemake.input.outdir 
    preds = predict_taxa(input_df, models)
    
    # output_file_path = os.path.join(out_dir, filename)
    preds.to_csv(os.path.join(out_dir, 'taxa_predictions.csv'), index=False)
    dists = get_distances(input_df, p_a_matrix)
    dists.to_csv(os.path.join(out_dir,'eucl_distances.csv'))
