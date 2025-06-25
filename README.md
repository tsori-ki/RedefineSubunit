# Protein Subunit Redefinition Pipeline

A computational pipeline for redefining protein subunits using structural confidence metrics and assembling complexes through graph-based approaches.

## Key Workflows
1. **AlphaFold Processing**: 
   - Multiple AlphaFold runs with confidence score (pLDDT) analysis
   - Subunit segmentation based on confidence thresholds (pLDDT >40)
   
2. **Graph Construction** (`complex_graph.py`):
   - NetworkX-based graph creation
   - PAE matrix analysis for edge definition
   - Chain-aware residue mapping and renaming

3. **Subunit Merging** (`redefine.py`):
   - Graph merging across multiple predictions
   - Sequence validation and conflict detection
   - High/low confidence segment classification

4. **HPC Integration**:
   - Slurm batch job management
   - Cluster-ready preprocessing workflows

## Technologies
- **Core**: Python 3, AlphaFold, CombFold
- **Libraries**: NetworkX, Biopython, NumPy
- **HPC**: Linux, Slurm
- **Analysis**: TM-score validation, MSA processing

## Quick Start
```bash
# Install dependencies
pip install -r requirements.txt

# Process AlphaFold outputs
python redefine.py /path/to/af_results/ chain_mapping.json subunits_info.json
```

## Project Structure
```
RedefineSubunit/
├── complex_graph.py       # Graph construction from structural data
├── redefine.py            # Core merging/validation logic
├── filter_high_subunits.py # Confidence-based filtering
└── preprocess/            # MSA generation & data formatting
