# Probe Design Script

This script designs tiled primer/probe amplicons from a SnapGene `.dna` file. It is intended for generating ordered primer sequences across a donor, construct, or target region for LOCK-seq probe design.

## What it does

The script:

- Reads a SnapGene `.dna` file
- Extracts the full sequence and annotated features
- Designs forward and reverse primer pairs across the sequence
- Applies basic primer filters:
  - primer length 20–24 bp
  - GC content 30–80%
  - avoids 4+ bp homopolymers
  - requires primers to end in G or C
  - balances forward/reverse primer GC content
- Targets amplicons around ~226 bp
- Spaces amplicons at least 650 bp apart
- Adds universal tag sequences to primers for ordering
- Optionally masks repeats using RepeatMasker
- Optionally masks ambiguous `N` bases

## Requirements

Python packages:

```bash
pip install pandas biopython snapgene_reader

Usage:
python probe_design.py PROJECT_NAME input_file.dna

Option with RepeatMasker:
python probe_design.py MyConstruct donor_construct.dna --use_repeatmasker --species mouse

Run with N-base masking:
python probe_design.py MyConstruct donor_construct.dna --mask_ns

| Argument             | Description                                             |
| -------------------- | ------------------------------------------------------- |
| `project_name`       | Prefix used for output files and primer names           |
| `input_file`         | Input SnapGene `.dna` file                              |
| `--use_repeatmasker` | Optional. Runs RepeatMasker before primer design        |
| `--species`          | Species passed to RepeatMasker. Default is `mouse`      |
| `--mask_ns`          | Optional. Treats ambiguous `N` bases as masked sequence |


Output

The main output is:

PROJECT_NAME_probes.csv

This file includes:

Primer name
Full primer sequence to order
Amplicon size
Amplicon GC%
Universal tag sequence
Forward or reverse primer sequence
Primer GC%
Primer length

If RepeatMasker is used, the script also writes:

PROJECT_NAME_masked_sequence.fa


| Column            | Description                                    |
| ----------------- | ---------------------------------------------- |
| `Primer Name`     | Primer ID using the project name and F/R index |
| `Primer to Order` | Universal tag + primer sequence                |
| `Placeholder`     | Empty column for downstream ordering formats   |
| `Amplicon Size`   | Size of the predicted amplicon                 |
| `Amplicon GC%`    | GC content of the amplicon                     |
| `Flag`            | Universal tag sequence                         |
| `Forward Primer`  | Forward primer sequence                        |
| `Reverse Primer`  | Reverse primer sequence                        |
| `Primer GC%`      | GC content of primer                           |
| `Primer Length`   | Primer length in bp                            |

Notes
Input must be a SnapGene .dna file.
RepeatMasker must be installed and available in the command line path if --use_repeatmasker is used.
The script currently uses fixed design parameters for amplicon size, primer length, spacing, and GC filtering.
Universal adapter/tag sequences are hard-coded in the script.
Features annotated in the SnapGene file are parsed and printed during runtime.

