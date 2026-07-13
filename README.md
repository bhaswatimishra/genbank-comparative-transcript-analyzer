# GenBank Comparative Transcript Analyzer

A Python-based bioinformatics pipeline for retrieving transcript records from NCBI GenBank, comparing transcript and protein characteristics across multiple species, and generating a phylogenetic tree using sequence alignment.

---

## Features

* Retrieve MANE Select transcript records from NCBI GenBank
* Convert GenBank records to FASTA format
* Extract transcript information:

  * Sequence length
  * GC percentage
  * Chromosome number
  * Topology
  * Exon count
* Analyze translated protein properties:

  * Molecular weight
  * Protein length
  * Isoelectric point (pI)
  * Aromaticity
* Compare multiple species in a CSV report
* Generate a multi-sequence FASTA file
* Perform multiple sequence alignment using MUSCLE
* Construct a phylogenetic tree using FastTree
* Display the resulting tree in ASCII format using Biopython

---

## Workflow

1. Enter the gene name.
2. The script automatically retrieves the human MANE Select transcript from NCBI.
3. A FASTA file is generated for the reference transcript.
4. Run BLAST manually using the generated FASTA sequence.
5. Enter accession IDs of homologous transcripts from other species.
6. The script downloads GenBank records for each selected species.
7. Transcript and protein properties are extracted.
8. A comparative CSV report is generated.
9. Protein sequences are aligned using MUSCLE.
10. A phylogenetic tree is generated using FastTree.

---

## Technologies Used

* Python
* Biopython
* Pandas
* NCBI Entrez API
* MUSCLE
* FastTree

---

## Requirements

### Python packages

```bash
pip install biopython pandas
```

### External software

The following programs must be installed separately and available in your system PATH:

* MUSCLE (v5 or later)
* FastTree

---

## Output Files

The pipeline generates:

* `transcript_comparator.csv` — Comparative transcript and protein analysis
* `all_sample.fasta` — Protein sequences of all selected species
* `aligned.fasta` — Multiple sequence alignment
* `tree.nwk` — Phylogenetic tree in Newick format

---

## Current Limitations

* BLAST search is **not automated**. Users must perform BLAST manually and provide accession IDs.
* The script currently uses **MANE Select** transcripts only for the human reference sequence.
* Limited input validation and error handling.
* Sequence alignment and tree construction require external software (MUSCLE and FastTree) to be installed.
* Tree visualization is currently limited to ASCII output in the terminal.
---

## Repository Structure

```
├── genbank_analyser.py
├── README.md
├── requirements.txt
├── transcript_comparator.csv (example output)
├── aligned.fasta (example output)
├── tree.nwk (example output)
```

---

This project was developed as a learning project to explore comparative genomics, sequence analysis, and phylogenetic analysis using Python and publicly available biological databases.
