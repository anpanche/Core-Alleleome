import logging
import os
import shutil
from pathlib import Path

import pandas as pd


def process_core_genes(
    pangenome_alignments_dir_path, alleleome_dir_path, pangene_summary_csv=None
):
    """
    Processes core genes by copying gene sequence files and creating output directories for each gene.

    Parameters:
    - main_dir: Path to the main directory containing gene folders.
    - alleleome_path: Path to the directory containing the gene data CSV files.
    """
    try:
        logging.info("Starting process_core_genes in QCQA_4")

        # Read core genes data
        alleleome_dir_path = Path(alleleome_dir_path)
        alleleome_dir_path.mkdir(parents=True, exist_ok=True)

        if pangene_summary_csv is None:
            pangene_summary_csv = alleleome_dir_path / "df_pangene_summary_v2.csv"
        else:
            pangene_summary_csv = Path(pangene_summary_csv)

        assert (
            pangene_summary_csv.is_file()
        ), f"Cannot find pangene_summary table at {pangene_summary_csv}"

        df = pd.read_csv(pangene_summary_csv)

        core_gene_list = df[df["pangenome_class_2"] == "Core"]["Gene"].unique().tolist()

        # Read core gene alleles data
        df1 = pd.read_csv(
            os.path.join(
                alleleome_dir_path,
                "core_alleles_with_length_less_than_2std_less_than_mean_length.csv",
            )
        )
        edit_list = df1["Gene"].tolist()
        new_gene_list = set(edit_list)

        # Process each core gene
        for gene in core_gene_list:
            # Copy gene sequence files if gene not in new_gene_list
            if gene not in new_gene_list:
                aa_allele_path = pangenome_alignments_dir_path + gene + "/input/"
                aa_file = aa_allele_path + "pangenes.faa"
                new_file = aa_allele_path + "pan_genes.faa"
                na_file = aa_allele_path + "pangenes.fna"
                new_na_file = aa_allele_path + "pan_genes.fna"
                shutil.copyfile(aa_file, new_file)
                shutil.copyfile(na_file, new_na_file)

            # Create output directory for the gene
            output_path = os.path.join(pangenome_alignments_dir_path, gene, "output")
            os.makedirs(output_path, exist_ok=True)
        logging.info("Completed process_core_genes in QCQA_4")
    except Exception as e:
        logging.error(f"Error in process_core_genes in QCQA_4: {e}")
        raise
