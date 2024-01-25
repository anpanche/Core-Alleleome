import logging
import subprocess
from pathlib import Path

import pandas as pd


def nucleotide_seq_align(
    pangenome_alignments_dir_path,
    alleleome_dir_path,
    pangene_summary_csv=None,
    blast_path="./resources/ncbi-blast-2.14.0+/bin/blastn",
):
    try:
        logging.info("Starting nucleotide_seq_align in nucleotide_sequence_alignment")
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

        core_na_query_list = (
            df["pangenome_class_2"].eq("Core").groupby(df["Gene"]).any()
        ).pipe(lambda x: x.index[x].tolist())

        if Path(blast_path).is_file():
            pass
        else:
            logging.warning(
                f"Cannot find blastn at {blast_path}. Trying default path..."
            )
            blast_path = "blastn"

        for r in range(len(core_na_query_list)):
            query = core_na_query_list[r]
            if pd.isnull(query) is False:
                out_file_name = (
                    pangenome_alignments_dir_path
                    / query
                    / "output"
                    / ("nucleotide_" + "blast_out_" + query + ".xml")
                )
                args = (
                    blast_path,
                    "-query",
                    pangenome_alignments_dir_path / query / "input" / "pan_genes.fna",
                    "-subject",
                    pangenome_alignments_dir_path
                    / query
                    / "input"
                    / ("nucleotide_consensus_" + query + ".fna"),
                    "-outfmt",
                    "5",
                )
                with open(out_file_name, "w+") as outfile:
                    subprocess.run(args, stdout=outfile)
        logging.info("Completed nucleotide_seq_align in nucleotide_sequence_alignment")
    except Exception as e:
        logging.error(
            f"Error in nucleotide_seq_align in nucleotide_sequence_alignment: {e}"
        )
        raise
