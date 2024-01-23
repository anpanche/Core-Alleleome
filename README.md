

# Alleleome: Core Gene Alleleome Generation for Pan-Genomes of Species

## Introduction
"Alleleome" is a specialized package designed to explore and analyze natural sequence variations within the Open Reading Frames (ORFs) of alleles of core genes in a species' pan-genome, both at the amino acid and nucleotide levels. This first version focuses on analyzing only the core genes of a pan-genome. It identifies variants such as substitutions, insertions, and deletions through a series of steps:

1. Initial QCQA of sequences.
2. Building consensus for each gene's allele set.
3. Pairwise alignment of consensus sequences with individual alleles.
4. Identification and generation of amino acid variant datasets.
5. Analysis of synonymous and non-synonymous substitutions from codons and corresponding amino acid data.

The Alleleome workflow is specifically tailored to study the natural sequence variations in core genes of the pan-genome of a species, with an emphasis on variations at the amino acid and nucleotide level.

### Publication
For more detailed information, refer to our publication:
[Early Release on BioRxiv](https://www.biorxiv.org/content/biorxiv/early/2023/09/22/2023.09.22.558971.full.pdf)

## Table of Contents
- [Introduction](#introduction)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Features](#features)
- [Built With](#built-with)
- [Versioning](#versioning)
- [Authors and Acknowledgment](#authors-and-acknowledgment)
- [License](#license)

## Getting Started

### Prerequisites
- Alleleome is tested and confirmed for Linux systems with the Conda package manager.
- Requires Python version 3.10 or higher.
- For optimal performance, especially when processing a large dataset, such as 1400 core genes and their respective alleles across 3400 strains, a high RAM capacity is strongly recommended.

### Installation
#### Using GitHub Repository
1. Clone the repository:
   ```
   git clone https://github.com/anpanche/PAN-Alleleome.git
   ```
2. Navigate to the Alleleome directory:
   ```
   cd PAN-Alleleome
   ```
3. Install the package:
   ```
   pip install .
   ```

## Usage

### Test Sample Data Run
Run Alleleome with sample data using:
   ```
   Alleleome
   ```
Refer to the `sample_data` directory for output organization and details.

### Running Your Species Core Genes
To analyze your species data:
   ```
   Alleleome --path1 path/to/pangenome_alignments --path2 path/to/alleleome
   ```


## Features
Alleleome introduces the concept of "ORF alleleome," encapsulating the gene alleles found across all strains of a species, thus providing a comprehensive view of genome-scale sequence variations. This analysis can be instrumental in understanding sequence diversity characteristics and natural selection processes across different species within a family. The study of the alleleome offers insights into the genetic basis of natural selection in a species.

Key features include:
- Analysis of sequence variants using the consensus sequence of ORFs.
- Identification of dominant amino acids and their variants at specific positions.
- Revealing natural sequence and structural variations compared with the consensus sequence and structural attributes.
- Identification of genome-scale synonymous and non-synonymous mutations through the analysis of codon changes and their corresponding amino acid changes."

## Built With
- Python and Biopython.
- Integrated with "BGCflow" workflow using SnakeMake.

## Versioning
- [Versioning details here]

## Authors and Acknowledgment
- [List of authors and contributors]
- Special thanks to [acknowledgments].

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

