from . import data_processing
from . import consensus_generate
from . import amino_acid_seq_alignment
from . import nucleotide_seq_alignment
from . import amino_acid_vars_generate
from . import codon_mutations
import argparse

def main():
    parser=argparse.ArgumentParser(description="Alleleome")
    parser.add_argument("--path1", help="Path argument 1")
    parser.add_argument("--path2", help="Path argument 2")
    args=parser.parse_args()
    if args.path1 and args.path2:
        data_processing.QC_QA(args.path1,args.path2)
        consensus_generate.Build_consensus(args.path1,args.path2)
        amino_acid_seq_alignment.Amino_acid_seq_align(args.path1,args.path2)
        amino_acid_vars_generate.Generate_amino_acid_vars(args.path1,args.path2)
        nucleotide_seq_alignment.Nucleotide_seq_align(args.path1,args.path2)
        codon_mutations.Codon_mut(args.path1,args.path2)
    else:
        args.path1='./Alleleome/sample_data/Oenococcus_oeni/pangenome_alignments/'
        args.path2='./Alleleome/sample_data/Oenococcus_oeni/alleleome/'
        data_processing.QC_QA(args.path1,args.path2)
        consensus_generate.Build_consensus(args.path1,args.path2)
        amino_acid_seq_alignment.Amino_acid_seq_align(args.path1,args.path2)
        amino_acid_vars_generate.Generate_amino_acid_vars(args.path1,args.path2)
        nucleotide_seq_alignment.Nucleotide_seq_align(args.path1,args.path2)
        codon_mutations.Codon_mut(args.path1,args.path2)