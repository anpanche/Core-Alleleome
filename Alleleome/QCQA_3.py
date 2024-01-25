import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def process_sequences(pangenome_alignments_dir_path, alleleome_dir_path):
    """
    Processes both nucleotide and amino acid sequence files based on specific criteria.

    Parameters:
    pangenome_alignments_dir_path: Path to the pangenome alignments directory containing gene folders.
    alleleome_dir_path: Path to the alleleome directory containing the CSV file with gene data required for further processing.
    """
    try:
        logging.info("Starting process_sequences in QCQA_3")
        pangenome_alignments_dir_path = Path(pangenome_alignments_dir_path)
        alleleome_dir_path = Path(alleleome_dir_path)
        df = pd.read_csv(
            alleleome_dir_path
            / "core_alleles_with_length_less_than_2std_less_than_mean_length.csv"
        )
        gene_locus_list = set(df["Gene"].to_list())
        locus_list = df["Locus_tag"].to_list()

        for file_type in ["fna", "faa"]:
            for gene in gene_locus_list:
                allele_path = pangenome_alignments_dir_path / gene / "input"
                seq_file = allele_path / f"pangenes.{file_type}"
                new_file = allele_path / f"pan_genes.{file_type}"
                with open(seq_file, "r") as file_in, open(new_file, "w") as file_out:
                    for seq_record in SeqIO.parse(file_in, "fasta"):
                        if seq_record.id not in locus_list:
                            desc = seq_record.description
                            allele_seq = str(seq_record.seq)
                            alleles = f">{desc}\n{allele_seq}\n"
                            file_out.write(alleles)
        logging.info("Completed process_sequences in QCQA_3")
    except Exception as e:
        logging.error(f"Error in process_sequences in QCQA_3: {e}")
        raise
