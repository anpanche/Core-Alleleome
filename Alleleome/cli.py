import argparse
import datetime
import logging
import os

from . import (
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

log_directory = "./log"
os.makedirs(log_directory, exist_ok=True)

current_time = datetime.datetime.now()
log_filename = os.path.join(
    log_directory, f"alleleome_{current_time.strftime('%Y-%m-%d_%H:%M:%S')}.log"
)

logging.basicConfig(
    filename=log_filename,
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def main():
    logging.info("Application started")
    try:
        parser = argparse.ArgumentParser(description="Alleleome")
        parser.add_argument("--path1", help="Path argument 1")
        parser.add_argument("--path2", help="Path argument 2")
        parser.add_argument(
            "--table", help="Path to custom CSV pangene summary table", default=None
        )
        parser.add_argument(
            "--log_to_terminal",
            action="store_true",
            help="Log output to terminal instead of file",
        )
        args = parser.parse_args()

        if args.log_to_terminal:
            for handler in logging.root.handlers[:]:
                logging.root.removeHandler(handler)
            logging.basicConfig(
                level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
            )
        else:
            pass

        if args.path1 and args.path2:
            QCQA_1.process_nucleotide_sequences(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            QCQA_2.analyze_gene_lengths(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            QCQA_3.process_sequences(args.path1, args.path2)
            QCQA_4.process_core_genes(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            build_consensus_sequence.build_consensus(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            amino_acid_sequence_alignment.amino_acid_seq_align(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            generate_amino_acid_variants.generate_amino_acid_vars(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            nucleotide_sequence_alignment.nucleotide_seq_align(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            codon_mutations.codon_mut(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
        else:
            args.path1 = "./Alleleome/sample_data/Oenococcus_oeni/pangenome_alignments/"
            args.path2 = "./Alleleome/sample_data/Oenococcus_oeni/alleleome/"
            QCQA_1.process_nucleotide_sequences(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            QCQA_2.analyze_gene_lengths(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            QCQA_3.process_sequences(args.path1, args.path2)
            QCQA_4.process_core_genes(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            build_consensus_sequence.build_consensus(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            amino_acid_sequence_alignment.amino_acid_seq_align(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            generate_amino_acid_variants.generate_amino_acid_vars(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            nucleotide_sequence_alignment.nucleotide_seq_align(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
            codon_mutations.codon_mut(
                args.path1, args.path2, pangene_summary_csv=args.table
            )
        logging.info("Application finished successfully")
    except Exception as e:
        logging.error(f"Application encountered an error: {e}", exc_info=True)
