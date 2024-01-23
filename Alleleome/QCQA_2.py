import os
from Bio import SeqIO
import pandas as pd
import logging

def analyze_gene_lengths(pangenome_alignments_dir_path,alleleome_dir_path):
    try:
        logging.info("Starting analyze_gene_lengths in QCQA_2")

        # Read gene summary
        df = pd.read_csv(os.path.join(alleleome_dir_path, 'df_pangene_summary_v2.csv'))

        # Extract core genes
        core_gene_list = (df['pangenome_class_2'].eq('Core').groupby(df['Gene']).any()).pipe(lambda x: x.index[x].tolist())

        nuc_data = []

        for k in core_gene_list:
            s = 0
            gene_length = []
            i = 0
            nuc_allele_file = os.path.join(pangenome_alignments_dir_path, k, 'input')
            nuc_file = os.path.join(nuc_allele_file, 'pangenes.fna')

            with open(nuc_file, 'r+') as f:
                for seq_record in SeqIO.parse(nuc_file, "fasta"):
                    g_id = seq_record.id
                    l_gene = len(seq_record)
                    gene_length.append(l_gene)
                    i = i + 1
                    s = s + l_gene
                    nuc_len = {
                        "Gene": k,
                        "Locus_tag": g_id,
                        "Length_of_allele": l_gene,
                    }
                    nuc_data.append(nuc_len)

        nuc_df = pd.DataFrame(nuc_data)
        nuc_df.to_csv(os.path.join(alleleome_dir_path, 'Genes_nuc_length_of_alleles.csv'))
        new_df = nuc_df.groupby(['Gene']).Length_of_allele.agg({'mean','std'})
        new_df['mean']=new_df['mean'].apply(lambda x:round(x,0))
        new_df['std']=new_df['std'].apply(lambda x:round(x,2))
        new_df.to_csv(alleleome_dir_path +'core_nuc_alleles_with_mean_std.csv')
        df2 = nuc_df.merge(new_df, left_on=['Gene'], right_index=True)
        df2.to_csv(os.path.join(alleleome_dir_path, 'core_nuc_alleles_with_locus_mean_std.csv'))
        df2['mean_2std']=(df2['mean']- 2*df2['std'])
        df3=df2[df2.Length_of_allele < df2.mean_2std]
        df3.to_csv(os.path.join(alleleome_dir_path,'core_alleles_with_length_less_than_2std_less_than_mean_length.csv'))
        logging.info("Completed analyze_gene_lengths in QCQA_2")   
    except Exception as e:
        logging.error(f"Error in analyze_gene_lengths in QCQA_2: {e}")         
        raise
