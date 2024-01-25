# Alleleome generation step II
# Alignment of each allele sequnces (amino acid sequence) of each gene with consensus sequence using BLASTp.
import logging
import subprocess
from pathlib import Path

import pandas as pd


def amino_acid_seq_align(
    pangenome_alignments_dir_path,
    alleleome_dir_path,
    pangene_summary_csv=None,
    blast_path="./resources/ncbi-blast-2.14.0+/bin/blastp",
):
    try:
        logging.info("Starting amino_acid_seq_align in amino_acid_sequence_alignment")
        pangenome_alignments_dir_path = Path(pangenome_alignments_dir_path)
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

        core_aa_query_list = (
            df["pangenome_class_2"].eq("Core").groupby(df["Gene"]).any()
        ).pipe(lambda x: x.index[x].tolist())

        if Path(blast_path).is_file():
            pass
        else:
            logging.warning(
                f"Cannot find blastp at {blast_path}. Trying default path..."
            )
            blast_path = "blastp"

        for r in range(len(core_aa_query_list)):
            query = core_aa_query_list[r]
            if pd.isnull(query) is False:
                out_file_name = (
                    pangenome_alignments_dir_path
                    / query
                    / "output"
                    / ("amino_acid_" + "blast_out_" + query + ".xml")
                )
                args = (
                    blast_path,
                    "-query",
                    pangenome_alignments_dir_path / query / "input" / "pan_genes.faa",
                    "-subject",
                    pangenome_alignments_dir_path
                    / query
                    / "input"
                    / ("amino_acid_consensus_" + query + ".faa"),
                    "-outfmt",
                    "5",
                )
                with open(out_file_name, "w+") as outfile:
                    subprocess.run(args, stdout=outfile)
        logging.info("Completed  amino_acid_seq_align in amino_acid_sequence_alignment")
    except Exception as e:
        logging.error(
            f"Error in amino_acid_seq_align in amino_acid_sequence_alignment: {e}"
        )
        raise
