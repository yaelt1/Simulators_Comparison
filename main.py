import json
import logging
from compare_simulators import compare_main
import os
import pandas as pd
from globals import feature_dict
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse

def file_to_array(paths:list[str], n:int=10000)-> np.array:
    data = np.zeros(len(paths)*n)
    i=0
    for path in paths:
        with open(path, "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                data[i] = (float(line))
                i+=1
    return data

def df_of_all(seq_length, indel_rate, result_path, modules, number_of_simulations):
    cols = [v for k,v in feature_dict.items()]
    cols.sort()
    cols.append("Module")
    feature_df = pd.DataFrame(columns= cols)
    for feature_number, feature_name in feature_dict.items(): 
        files=[]
        for mod in modules:
            files.append(os.path.join(result_path,mod.upper(),feature_name+".txt"))
        arr = file_to_array(files, number_of_simulations)
        feature_df[f"{feature_name}"] = arr
    
    samples = 0
    for mod in modules:        
        feature_df["Module"][samples:samples+number_of_simulations] = mod
        samples += number_of_simulations
    
    annotation_str = f"Indel Rate:{indel_rate}\nSequence Length:{seq_length}"
    for feature_number, feature_name  in feature_dict.items(): 
        sns.histplot(data=feature_df, x=f'{feature_name}', hue='Module', alpha=0.7)
        plt.title(f"{feature_name}")
        plt.annotate(annotation_str, (1, 2), xytext=(0.86, 1.005), fontsize=9, ha='center', xycoords='axes fraction')
        fig_path = os.path.join(result_path, f"{feature_name}_hist.png")
        plt.savefig(fig_path)
        plt.show()
def parse_args():
    parser = argparse.ArgumentParser(description='Simulator Comparison Tool')
    parser.add_argument('--config', type=str, help='Path to the configuration JSON file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if args.config:
        with open(args.config) as f:
            data = json.load(f)
    else:
        print("Error: Please provide a configuration JSON file using --config option.")
        exit(1)
    tree_filepath = data['tree_filepath']
    num_nodes = data['num_nodes']
    indel_rate = data['indel_rate']
    length_seq = data['length_seq']
    result_path =  data['result_path']
    modules = data['modules']
    number_of_simulations = data['number_of_simulations']
    features = set(data['features'])
    result_path = os.path.join(result_path, f"{num_nodes}Nodes_{indel_rate}Rate_{length_seq}Length")
    os.makedirs(result_path, exist_ok=True)
    log_path = os.path.join(result_path, 'Simulators_comparison.log')
    logging.basicConfig(filename=log_path, level=logging.INFO) 
    logging.info('Parameters read from the json file')
    logging.info('tree_filepath: %s', tree_filepath)
    logging.info('num_nodes: %d', num_nodes)
    logging.info('indel_rate: %f', indel_rate)
    logging.info('length_seq: %d', length_seq)
    logging.info('result_path: %s', result_path)
    logging.info('Starting the comparison of simulators')
    compare_main(indel_rate, length_seq, result_path, tree_filepath, num_nodes, features, modules, number_of_simulations)    
    df_of_all(length_seq, indel_rate, result_path, modules, number_of_simulations)
    
        
        
