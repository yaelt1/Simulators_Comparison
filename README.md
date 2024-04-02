The Simulator Comparison Tool is a project designed to compare different simulation tools for sequence generation and analysis.

## Installation
1. Clone the repository.
3. Set up the project configuration by editing `simulator.json"

The `simulator.json` file contains configuration parameters for the project. Example content:
```json
{
    "tree_filepath": "/path/to/tree/file",
    "num_nodes": 100,
    "indel_rate": 0.05,
    "length_seq": 500,
    "result_path": "/path/to/results",
    "modules": ["alisim", "indelible", "sailfish"],
    "number_of_simulations": 1000,
    "features": [1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
}

## Global Variables

### Indelible Configuration
- `INDELIBLE_EXECUTABLE`: Path to the INDELible executable file.

### Feature Dictionary
- `feature_dict`: Dictionary mapping feature numbers to their corresponding names.
```python
feature_dict = {
    1: "Number of Gaps",
    2: "Avg Gap Length",
    5: "Length 1 Gaps",
    6: "Length 2 Gaps",
    7: "Length 3 Gaps",
    8: "Length 4+ Gaps",
    9: "Alignment Length",
    10: "Shortest Sequences",
    11: "Longest Sequences",
    12: "Cols Without Gaps",
    13: "Cols with 1 Gap",
    14: "Cols With 2 Gaps",
    15: "Cols With n-1 Gaps"
}
