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
from create_figs import craete_subfigs
      
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
    # config = "simulator_real.json"
    # with open(config) as f:
    #         data = json.load(f)
    tree_filepath = data['tree_filepath']
    num_nodes = data['num_nodes']
    indel_rate = data['indel_rate']
    length_seq = data['length_seq']
    result_path =  data['result_path']
    modules = data['modules']
    number_of_simulations = data['number_of_simulations']
    features = set(data['features'])
    name = f"{num_nodes}Nodes_{indel_rate}Rate_{length_seq}Length"
    result_path = os.path.join(result_path, name)
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
    craete_subfigs(length_seq, indel_rate, result_path, modules, number_of_simulations,name )
        
        
