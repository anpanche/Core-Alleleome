#Alleleome generation step II- Alignment of each allele sequnces (amino acid sequence) of each gene with consensus sequence using BLASTp.
import pandas as pd
import os
import subprocess


def Amino_acid_seq_align(pangenome_alignments_dir_path,alleleome_dir_path):
    df=pd.read_csv(alleleome_dir_path +  'df_pangene_summary_v2.csv')
    core_aa_query_list = (df['pangenome_class_2'].eq('Core').groupby(df['Gene']).any()).pipe(lambda x:x.index[x].tolist())
    blast_path='./resources/ncbi-blast-2.14.0+/bin/blastp'
    for r in range(len(core_aa_query_list)):
        query=core_aa_query_list[r]
        if (pd.isnull(query)==False):
            out_file_name = pangenome_alignments_dir_path + query + '/output/' + 'amino_acid_' + 'blast_out_' + query + '.xml'
            args = (blast_path,
                    '-query', pangenome_alignments_dir_path + query + '/input/'+'pan_genes.faa',
                    '-subject', pangenome_alignments_dir_path + query + '/input/'+ 'amino_acid_consensus_'+ query +'.faa' ,
                    '-outfmt', '5',
                )
            with open(out_file_name, 'w+') as outfile:
                subprocess.run(args, stdout=outfile)


