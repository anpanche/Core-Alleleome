import glob
from Bio.Blast import NCBIXML
import pandas as pd
from Bio.Seq import Seq
import os



table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        '---':'-'
        }



def Codon_mut(pangenome_alignments_dir_path,alleleome_dir_path):
    gene_syno_nonsyno_df = pd.DataFrame()
    df=pd.read_csv(alleleome_dir_path +  'df_pangene_summary_v2.csv')
    core_gene_list = (df['pangenome_class_2'].eq('Core').groupby(df['Gene']).any()).pipe(lambda x:x.index[x].tolist())
    for blast_out_file in core_gene_list:
        blast_output_file_path = pangenome_alignments_dir_path + blast_out_file + '/output/' 
        blast_output_file = blast_output_file_path +'nucleotide_blast_out_'+ blast_out_file + '.xml'
        gene = blast_out_file.replace('nucleotide_blast_out_','').replace('.xml','')
        for record in NCBIXML.parse(open(blast_output_file)):
            if len(record.alignments) > 0:
                #Description of available members: https://biopython.org/docs/1.75/api/Bio.Blast.Record.html
                subject_match_start_pos = record.alignments[0].hsps[0].sbjct_start
                blast_subject_str = record.alignments[0].hsps[0].sbjct
                blast_match_str = record.alignments[0].hsps[0].match
                blast_query_str = record.alignments[0].hsps[0].query
                ref_id= record.alignments[0].hit_id
                ref_info= record.alignments[0].hit_def # new added#### column name to be added####
                ref_info=ref_info.split(" ", 1)
                ref_id=ref_info[0]
                query_info= record.query
                query_info=query_info.split("|", 1)
                query_description=query_info[0]
                query_description=query_description.split(" ", 1)
                query_id=query_description[0]
                query_desc=query_description[1]
                GCF_id=query_info[1]
                blast_query_str=blast_query_str.replace('-','N')
                blast_query_str=Seq(blast_query_str)
                blast_subject_str=blast_subject_str.replace('-','N')
                blast_subject_str=Seq(blast_subject_str)
                end_sub = len(blast_subject_str) - (len(blast_subject_str) %3) - 1
                end_query= len(blast_query_str) - (len(blast_query_str) %3) - 1
                for j in range(0,end_sub,3):
                    codon_s = blast_subject_str[j:j+3]
                    codon_q = blast_query_str[j:j+3]
                    if codon_s == codon_q:
                        pass
                    elif codon_s in table:
                        aa_s = codon_s.translate()
                        aa_q = codon_q.translate()
                        if aa_s == aa_q:
                            df = pd.DataFrame.from_dict({
                                'Gene': gene,
                                'Cons_codon': codon_s,
                                'Query_codon': codon_q,
                                "Codon_position": j,
                                "Cons_aa": aa_s,
                                "Query_aa" : aa_q,
                                "AA_mutation_effect": "Synonymous",
                                "Cons_locus_tag": ref_id,
                                "Query_locus_tag" : query_id,
                                "Query_description" : query_desc,
                                "GCF_id" : GCF_id,
                            }, orient='index')  # Can't use columns due to diff parsing of 'AA range' input
                            df = df.T
                            gene_syno_nonsyno_df = pd.concat([df, gene_syno_nonsyno_df])
                        else:
                            df = pd.DataFrame.from_dict({
                                'Gene': gene,
                                'Cons_codon': codon_s,
                                'Query_codon': codon_q,
                                "Codon_position": j,
                                "Cons_aa": aa_s,
                                "Query_aa" : aa_q,
                                "AA_mutation_effect": "Non-synonymous",
                                "Cons_locus_tag": ref_id,
                                "Query_locus_tag" : query_id,
                                "Query_description" : query_desc,
                                "GCF_id" : GCF_id,
                            }, orient='index')  # Can't use columns due to diff parsing of 'AA range' input
                            df = df.T
                            gene_syno_nonsyno_df = pd.concat([df, gene_syno_nonsyno_df]) # Accounting for multiples of same allele
                    else:
                        aa_s = 'X'
                        aa_q = 'X'

        gene_syno_nonsyno_df.to_csv(os.path.join(alleleome_dir_path,'pan_core_gene_syno_non_syno_df.csv'))

