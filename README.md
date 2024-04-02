The Simulator Comparison Tool is a project designed to compare different simulation tools for sequence generation and analysis.

## Installation
1. Clone the repository.
3. Set up the project configuration by editing `simulator.json"

## Running the Program
To run the Simulator Comparison Tool, execute the following command from the project's main directory:
```python
python main.py --config /path/to/simulator.json
```
## Explanation of Parameters in Config file
- tree_filepath: Path to the file containing the tree structure used for simulation.
- num_nodes: The number of nodes in the tree.
- indel_rate: The rate of indel (insertion and deletion) events during sequence simulation.
- length_seq: The length of the sequences to be generated.
- result_path: Path to the directory where simulation results will be saved.
- modules: List of simulation modules to be used (e.g., "alisim", "indelible", "sailfish").
- number_of_simulations: The number of simulations to perform for each module.
- features: List of feature numbers to analyze in the simulation results (e.g., 1 for "Number of Gaps", 2 for "Avg Gap Length", etc.). Explanation of each feature number can be found in the Global Variables section.

`simulator.json` Example Content:
```json
{
    "tree_filepath": "/path/to/treefile",
    "num_nodes": number_of_nodes_in_tree,
    "indel_rate": 0.01,
    "length_seq": 500,
    "result_path": "/path/to/results",
    "modules": ["alisim", "indelible", "sailfish"],
    "number_of_simulations": 1000,
    "features": [1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
}
```

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
```


