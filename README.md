# RedefineSubunit

RedefineSubunit is a pipeline for predicting the structure of protein complexes by first identifying and redefining stable subunit boundaries using AlphaFold confidence scores. This research-oriented project focuses on data-driven subunit segmentation and reconstructing multimeric protein structures with improved interpretability and accuracy.

## Project Goal

The main objective is to build a better understanding of which parts of a protein are structurally stable — and to redefine subunits accordingly — before assembling the full protein complex. By leveraging confidence scores (pLDDT) from AlphaFold predictions, we can isolate robust fragments and use them to improve the accuracy and modularity of complex prediction.

## Motivation

Common protein-structure prediction algorithms often rely on predefined chain divisions. However, flexible or partially disordered regions can reduce prediction quality. Our approach dynamically redefines subunits based on how reliably they fold across AlphaFold runs, allowing for more biologically meaningful modeling and modularity.

## Key Features

- Analyze multiple AlphaFold outputs to extract high-confidence regions (pLDDT > 40).
- Automatically redefine protein subunits based on structural confidence.
- Combine both stable and low-confidence regions into candidate inputs for complex reconstruction.
- Assemble protein complexes from these redefined subunits via CombFold.
- Modular preprocessing steps for MSA generation, data cleaning, and postprocessing.

## Project Structure

```
RedefineSubunit/
│   cif_to_pdb.py
│   complex_graph.py
│   filter_high_subunits.py
│   redefine.py
│   README.md
│   requirements.txt
│   ...
└───preprocess/
        af3_json_input_to_fasta.py
        ... 
```

- Core logic is implemented in `redefine.py` and related scripts.
- Preprocessing handled in `preprocess/`, including MSA generation and data formatting.
- TM-score scripts are used for structural validation and evaluation.

## Technologies

- Python 3
- AlphaFold (Multimer model)
- Biopython
- NumPy
- External tools: MSA generation, TM-score

## Current Status

This is a research-stage project. Preliminary results show performance comparable to CombFold (without extended runs on subunit triplets/quartets), with expectations for improved accuracy upon full-scale testing.

## 📬 Contact

Feel free to reach out if you'd like to learn more about the project, collaborate, or discuss protein modeling challenges.

—
Developed as part of a research project in structural bioinformatics.