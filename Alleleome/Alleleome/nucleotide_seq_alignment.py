import pandas as pd
import os
import random
import pathlib
from os.path import join, getsize
import subprocess

def Nucleotide_seq_align(pangenome_alignments_dir_path,alleleome_dir_path):
    df=pd.read_csv(alleleome_dir_path + 'df_pangene_summary_v2.csv')
    core_na_query_list = (df['pangenome_class_2'].eq('Core').groupby(df['Gene']).any()).pipe(lambda x:x.index[x].tolist())
    #blast_path='/home/azureuser/datadrive/ncbi-blast-2.13.0+/bin/blastn'
    blast_path ='./resources/ncbi-blast-2.14.0+/bin/blastn'
    for r in range(len(core_na_query_list)):
        query=core_na_query_list[r]
        if (pd.isnull(query)==False):
            out_file_name = pangenome_alignments_dir_path + query + '/output/' + 'nucleotide_' + 'blast_out_' + query + '.xml'
            args = (blast_path,
                    '-query', pangenome_alignments_dir_path + query + '/input/'+'pan_genes.fna',
                    '-subject', pangenome_alignments_dir_path + query + '/input/'+ 'nucleotide_consensus_'+ query +'.fna' ,
                    '-outfmt', '5',
                )
            with open(out_file_name, 'w+') as outfile:
                subprocess.run(args, stdout=outfile)

