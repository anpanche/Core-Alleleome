##Data preprocessing step I- Calculate the presence of each core gene (nucleotide sequence) in the number of strains 
#Input:df_pangene_summary_v2.csv, pangenome_alignments (directory of genes- nucleotide sequence),
#output: Genes_nuc_length_number_of_strains.csv
import pkg_resources
import pandas as pd
import os
from Bio import SeqIO
import shutil

def process_file(input_file, output_file, locus_id_list):
    with open(input_file, 'r') as f, open(output_file, 'w') as f1:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            if seq_record.id not in locus_id_list:
                desc = seq_record.description
                allele_seq = str(seq_record.seq)
                alleles = ''.join(['>' + desc + '\n' + allele_seq + '\n'])
                f1.write(alleles)

def QC_QA(pangenome_alignments_dir_path,alleleome_dir_path):
    df=pd.read_csv(alleleome_dir_path +'df_pangene_summary_v2.csv')
    core_gene_list = (df['pangenome_class_2'].eq('Core').groupby(df['Gene']).any()).pipe(lambda x:x.index[x].tolist())
    genes_df = pd.DataFrame()
    for k in core_gene_list:
        s=0
        i=0
        nuc_allele_file = pangenome_alignments_dir_path + k + '/input/'
        nuc_file=  nuc_allele_file + 'pangenes.fna'
        f=open(nuc_file,'r+')
        for seq_record in SeqIO.parse(nuc_file, "fasta"):
            g_id=seq_record.id
            l_gene=len(seq_record)
            i=i+1
            s=s+l_gene
        df=pd.DataFrame.from_dict({
        "Gene":k,
        "Number_of_strains":i,
        }, orient='index')
        df=df.T
        genes_df=pd.concat([df,genes_df])
        genes_df.to_csv(os.path.join(alleleome_dir_path,'Genes_nuc_length_number_of_strains.csv'))

#Data preprocessing stepIII- Get the locus tag and length of each allele (nucleotide sequence) of the set of genes
#Input: df_pangene_summary_v2.csv, path to pangeome_alignments genes (.fna)
#Output: Genes_nuc_length_of_alleles.csv
#pangene_summary_filename=df_pangene_summary_v2.csv->This is an output file from pangenome analysis
    df=pd.read_csv(alleleome_dir_path +'df_pangene_summary_v2.csv')
    core_gene_list = (df['pangenome_class_2'].eq('Core').groupby(df['Gene']).any()).pipe(lambda x:x.index[x].tolist())
    genes_df = pd.DataFrame()
    for k in core_gene_list:
        s=0
        gene_length=[]
        i=0
        nuc_allele_file = pangenome_alignments_dir_path + k + '/input/'
        nuc_file=  nuc_allele_file + 'pangenes.fna'
        f=open(nuc_file,'r+')
        for seq_record in SeqIO.parse(nuc_file, "fasta"):
            g_id=seq_record.id
            l_gene=len(seq_record)
            gene_length.append(l_gene)
            i=i+1
            s=s+l_gene
            df=pd.DataFrame.from_dict({
                "Gene":k,
                "Locus_tag":g_id,    
                "Length_of_allele":l_gene,    
            }, orient='index')
            df=df.T
            genes_df=pd.concat([df,genes_df])
            genes_df.to_csv(os.path.join(alleleome_dir_path,'Genes_nuc_length_of_alleles.csv'))
#Data preprocessing stepIV- Calculate mean length, standard deviation of the length of the set of genes
#Input:Genes_nuc_length_of_alleles.csv
#Output:core_nuc_alleles_with_mean_std.csv,core_nuc_alleles_with_locus_mean_std.csv
    df1=pd.read_csv(alleleome_dir_path + 'Genes_nuc_length_of_alleles.csv')
    new_df1 = df1.groupby(['Gene']).Length_of_allele.agg({'mean','std'})
    new_df1['mean']=new_df1['mean'].apply(lambda x:round(x,0))
    new_df1['std']=new_df1['std'].apply(lambda x:round(x,2))
    new_df1.to_csv(os.path.join(alleleome_dir_path,'core_nuc_alleles_with_mean_std.csv'))
    df1 = df1.merge(new_df1, left_on=['Gene'], right_index=True)
    df1.to_csv(os.path.join(alleleome_dir_path,'core_nuc_alleles_with_locus_mean_std.csv'))
#Data preprocessing stepV-Get the genes with their alleles with length less than two standard deviation of the mean length 
#Input:core_nuc_alleles_with_locus_mean_std.csv,
#Output:core_alleles_with_length_less_than_2std_less_than_mean_length.csv
    df2=pd.read_csv(alleleome_dir_path +'core_nuc_alleles_with_locus_mean_std.csv')
    df2['mean_2std']=(df2['mean']- 2*df2['std'])
    df3=df2[df2.Length_of_allele < df2.mean_2std]
    df3.to_csv(os.path.join(alleleome_dir_path,'core_alleles_with_length_less_than_2std_less_than_mean_length.csv'))
#Data preprocessing stepVI-Removal of amino acid sequences with specific locus tag with length which is less than two standard deviation of the mean length from the pangenome and create a new file of final set of alleles (of amino acid sequences).In this step
#the new file is created only of those genes in which alleles with prescribed length conditions are matched. 
#Input: core_alleles_with_length_less_than_2std_less_than_mean_length.csv
    df4=pd.read_csv(alleleome_dir_path + 'core_alleles_with_length_less_than_2std_less_than_mean_length.csv')
    gene_locus_list=df4['Gene'].to_list()
    gene_locus_list=set(gene_locus_list)
    locus_list=df4['Locus_tag'].to_list()

    for s in gene_locus_list:
        if s in gene_locus_list:
            aa_allele_path = pangenome_alignments_dir_path + s + '/input/'
            aa_file = aa_allele_path + 'pangenes.faa'
            na_file = aa_allele_path + 'pangenes.fna'
            new_file = ''.join(aa_allele_path + 'pan_genes.faa')
            new_na_file = ''.join(aa_allele_path + 'pan_genes.fna')
            process_file(aa_file, new_file, locus_list)
            process_file(na_file, new_na_file, locus_list)
       

    df5=pd.read_csv(alleleome_dir_path + 'df_pangene_summary_v2.csv')
    core_gene_list = (df5['pangenome_class_2'].eq('Core').groupby(df5['Gene']).any()).pipe(lambda x:x.index[x].tolist())
    #dataframe of core gene alleles with genes not satisfying the condition to be removed
    df6=pd.read_csv(alleleome_dir_path + 'core_alleles_with_length_less_than_2std_less_than_mean_length.csv')
    edit_list=df6['Gene'].to_list()
    new_gene_list=set(edit_list)
    for f in core_gene_list:
        if f not in new_gene_list:
            aa_allele_path = pangenome_alignments_dir_path + f + '/input/' 
            #aa_file=  aa_allele_path + 'pangenes.faa'
            aa_file=aa_allele_path +'pangenes.faa'
            #new_file= aa_allele_path + 'pan_genes.faa'
            new_file=aa_allele_path +'pan_genes.faa'
            #na_file = aa_allele_path + 'pangenes.fna'
            na_file=aa_allele_path +'pangenes.fna'
            #new_na_file = aa_allele_path + 'pan_genes.fna'
            new_na_file=aa_allele_path +'pan_genes.fna'
            #os.rename(*(os.path.join(aa_allele_path, fname) for fname in (aa_file, new_file)))
            #os.rename(*(os.path.join(aa_allele_path, fname) for fname in (na_file, new_na_file)))
            shutil.copyfile(aa_file, new_file)
            shutil.copyfile(na_file, new_na_file)

    #Data preprocessing stepVIII-Create new output directory to save the results of the further processing.
    j=0
    #query_list_1=[name for name in os.listdir(MAIN_DIR_1) if os.path.isdir(os.path.join(MAIN_DIR_1, name))]
    for i in core_gene_list:
        root_dir = os.path.join(pangenome_alignments_dir_path + i +'/')
        dir="output"
        path = os.path.join(root_dir, dir)
        j=j+1
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path)





