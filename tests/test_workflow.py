import pytest

from Alleleome import (
    QCQA_1,
    QCQA_2,
    QCQA_3,
    QCQA_4,
    amino_acid_sequence_alignment,
    build_consensus_sequence,
    codon_mutations,
    generate_amino_acid_variants,
    nucleotide_sequence_alignment,
)

path1 = "./Alleleome/sample_data/Oenococcus_oeni/pangenome_alignments/"
path2 = "./Alleleome/sample_data/Oenococcus_oeni/alleleome/"
table = None  # Replace with your pangene_summary_csv if any


def test_process_nucleotide_sequences():
    QCQA_1.process_nucleotide_sequences(path1, path2, pangene_summary_csv=table)


def test_analyze_gene_lengths():
    QCQA_2.analyze_gene_lengths(path1, path2, pangene_summary_csv=table)


def test_process_sequences():
    QCQA_3.process_sequences(path1, path2)


def test_process_core_genes():
    QCQA_4.process_core_genes(path1, path2, pangene_summary_csv=table)


def test_build_consensus():
    build_consensus_sequence.build_consensus(path1, path2, pangene_summary_csv=table)


def test_amino_acid_seq_align():
    amino_acid_sequence_alignment.amino_acid_seq_align(
        path1, path2, pangene_summary_csv=table
    )


def test_generate_amino_acid_vars():
    generate_amino_acid_variants.generate_amino_acid_vars(
        path1, path2, pangene_summary_csv=table
    )


def test_nucleotide_seq_align():
    nucleotide_sequence_alignment.nucleotide_seq_align(
        path1, path2, pangene_summary_csv=table
    )


def test_codon_mut():
    codon_mutations.codon_mut(path1, path2, pangene_summary_csv=table)
