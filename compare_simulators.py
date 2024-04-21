import argparse
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os 
from Bio import Align, AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from alisim_runner import AlisimCommandline
from sailfish_gen import Sailfish_sim
from indelible_dict import indelible_w_indels
from indelible_runner import indelible_wrapper
from Bio.Align import MultipleSeqAlignment
import globals
import logging
from globals import feature_dict
class features_calculator: 
    def __init__(self, module) -> None:
        self.module = module
        self.num_of_gaps = 0
        self.total_length_of_gaps=0
        self.gaps_of_len_1 =0
        self.gaps_of_len_2 =0
        self.gaps_of_len_3 =0
        self.gaps_of_len_4p =0
        self.min_val = 0
        self.max_val = 0
        self.len_of_align = 0
        self.cols_without_gaps = 0
        self.cols_w_1_gap = 0
        self.cols_w_2_gaps = 0
        self.cols_w_n_minus_1 = 0
        

def colums_gap_calc(path, features: features_calculator):
    with open(path) as f:
        first_line = f.readline()
        file_format = "fasta" if ">" in first_line else "phylip"
    
    
    records = list(SeqIO.parse(path, file_format))
    n = len(records)
    alignment = Align.MultipleSeqAlignment(records)
    for col in range(alignment.get_alignment_length()):
        column_data = str(alignment[:, col])
        gap_counter = column_data.count("-")
        
        if gap_counter ==0:
            features.cols_without_gaps +=1
        if gap_counter == 1:
            features.cols_w_1_gap +=1
        if gap_counter ==2:
            features.cols_w_2_gaps +=1
        if gap_counter == n-1:
            features.cols_w_n_minus_1 +=1    
           

def gap_blocks_features_calc(path:str, features: features_calculator )->None:
    with open(path) as f:
        first_line = f.readline()
        file_format = "fasta" if ">" in first_line else "phylip"
        f.seek(0)
        
        for record in SeqIO.parse(f, file_format):
            features.len_of_align = len(record.seq)
            cur_gap = 0
            i=0
            while i <len(record.seq):
                if record.seq[i] == "-" :
                    cur_gap +=1
        
                elif cur_gap != 0 :
                    features.total_length_of_gaps += cur_gap
                    features.num_of_gaps += 1
                    if cur_gap >=4:
                            features.gaps_of_len_4p +=1
                    if cur_gap==1:
                            features.gaps_of_len_1 +=1
                    if cur_gap==2:
                            features.gaps_of_len_2 +=1
                    if cur_gap==3:
                            features.gaps_of_len_3 +=1
                    cur_gap = 0
                i+=1
                        

def longest_n_shortest(path:str, features: features_calculator)->tuple:
    """ 
    find longest and shortest length sequnces from simulator's result

    """
    min_val = math.inf
    max_val = - (math.inf)
    with open(path) as f:
        first_line = f.readline()
        file_format = "fasta" if ">" in first_line else "phylip"
        f.seek(0)
        for record in SeqIO.parse(f, file_format):
            if len(record.seq) < min_val:
                min_val = len(record.seq)
            if len(record.seq) > max_val:
                max_val = len(record.seq)
    features.min_val = min_val
    features.max_val= max_val


def get_len_without_gaps(seq:str)->int:
    cnt = 0
    cnt = sum([1 if seq[i]!= "-" else 0 for i in range(len(seq))])
    return cnt


def longest_n_shortest_v2(path:str, features: features_calculator)->tuple:
    """ 
    find longest and shortest length sequnces from simulator's result - output with gaps

    """
    min_val = math.inf
    max_val = - (math.inf)
    with open(path) as f:
        first_line = f.readline()
        file_format = "fasta" if ">" in first_line else "phylip"
        f.seek(0)
        for record in SeqIO.parse(f, file_format):
            len_without_gaps = get_len_without_gaps(record.seq)
            if len_without_gaps < min_val:
                min_val = len_without_gaps
            if len_without_gaps > max_val:
                max_val = len_without_gaps
    features.min_val = min_val
    features.max_val= max_val


def remove_all_gaps_columns(alignment):
    gap_columns = [i for i in range(alignment.get_alignment_length()) if all(record.seq[i] == "-" for record in alignment)]

    # Create a new alignment excluding columns with only gaps
    filtered_alignment = MultipleSeqAlignment(
    SeqRecord(Seq("".join(record.seq[i] for i in range(alignment.get_alignment_length()) if i not in gap_columns)), id=record.id)
    for record in alignment)
    return filtered_alignment


def write_filtered_fasta(filtered_alignment, output_filename):
    with open(output_filename, "w") as output_handle:
        AlignIO.write(filtered_alignment, output_handle, "fasta")    


def analyze_simulation_output(module:str , output_filename_with_gap:str=None, output_without_gap:str=None, features:list= [1,2,5,6,7,8,9,10,11,12,13,14,15])-> dict:
    """_summary_

    Args:
        module (str): module name
        output_filename_with_gap (_type_): output with gaps 
        output_without_gap (str, optional): output without gaps  - Defaults to None.
        features: featues numbers to check.

    Returns:
        dict: _description_
    """
    result_dict = {}
    # Create a feature Calculator object to keep the data
    features_holder = features_calculator(module)
    
    if module== "alisim":
        alignment = AlignIO.read(output_filename_with_gap, "fasta")
        
    elif module == "indelible":
        alignment = AlignIO.read(output_filename_with_gap, "phylip-relaxed")
    
    if module != "sailfish":
        filtered_alignment = remove_all_gaps_columns(alignment)
        write_filtered_fasta(filtered_alignment, output_filename_with_gap)
        # logging.info("Filtered alignment written to %s", output_filename_with_gap)
           
    # check if there are blocks-related features
    if features&{1,2,5,6,7,8,9}:
        gap_blocks_features_calc(output_filename_with_gap, features_holder)
        result_dict[1] = features_holder.num_of_gaps
        if features_holder.num_of_gaps != 0:
            result_dict[2] = features_holder.total_length_of_gaps / features_holder.num_of_gaps
        result_dict[5] = features_holder.gaps_of_len_1
        result_dict[6] = features_holder.gaps_of_len_2
        result_dict[7] = features_holder.gaps_of_len_3
        result_dict[8] = features_holder.gaps_of_len_4p
        result_dict[9] = features_holder.len_of_align
    
    # check if there are longest and shortest related features
    if features&{10,11}:
        if output_without_gap == None:
            longest_n_shortest_v2(output_filename_with_gap, features_holder)
        else:    
            longest_n_shortest(output_without_gap, features_holder)
        result_dict[10] = features_holder.min_val
        result_dict[11] = features_holder.max_val
    
    # check if there are columns-related features    
    if features&{12,13, 14, 15}:
        colums_gap_calc(output_filename_with_gap, features_holder)
        result_dict[12] = features_holder.cols_without_gaps
        result_dict[13]= features_holder.cols_w_1_gap
        result_dict[14] = features_holder.cols_w_2_gaps
        result_dict[15] = features_holder.cols_w_n_minus_1
    return result_dict


def get_features(module, output_path_with_gaps, output_path_without_gaps, features, all_results):
    features_result = analyze_simulation_output(module, output_path_with_gaps, output_path_without_gaps, features)
    for key in features_result.keys():
        if key in all_results:
            all_results[key].append(features_result[key])
        else:
            all_results[key] = [features_result[key]]
    return all_results
                    

def multi_time_simulation(tree_path:str, module:str, features:list,module_res_path:str ,seq_length:int=10000, indel_rate:float=0.01, number_of_simulations =1000)-> None:
    """
    Runs a module 10000 times, checks the features of interest and saves the results in txt file for each run
    
    :module - nodule name (AliSim, INDELible , Sailfish)
    
    :features - a list of features numbers needed to be check in each run

    """
    all_results = {}
    for i in range(number_of_simulations):
        if module.lower() == "indelible":
                indelible_wrapper(seq_length=seq_length ,indel_rate= indel_rate,result_path= module_res_path , indel_dict= indelible_w_indels, tree_path=tree_path)
                output_path_with_gaps =globals.indelible_with_gaps
                output_path_without_gaps = globals.indelible_without
        
        elif module.lower() == "alisim":
                AlisimCommandline(tree_path = tree_path,result_path= module_res_path ,seq_length= seq_length, indel_rate= indel_rate)
                output_path_with_gaps = globals.alisim_with_gaps
                output_path_without_gaps = globals.alisim_without
                
        elif module.lower() == "sailfish":
                os.chdir(module_res_path)
                Sailfish_sim(tree_path= tree_path,seq_length= seq_length, inse_r= indel_rate, del_r=indel_rate)
                output_path_with_gaps = globals.sailfish_with_gaps
                output_path_without_gaps = globals.sailfish_without
                
            
        else:        
            print("Model is not Supported")
            break
        all_results=get_features(module, output_path_with_gaps, output_path_without_gaps, features, all_results)
    return all_results


def save_results(all_results, path):
    cols = sorted(list(all_results.keys()))
    df = pd.DataFrame(columns=cols)
    for key in all_results.keys():
        df[key] = all_results[key]
    print(os.path.join(path, "features.csv"))
    df.to_csv(os.path.join(path, "features.csv"))
               

def compare_main(indel_rate:float, seq_length:int, result_path:str, tree_path:str, num_nodes:int, features:set, modules:list, number_of_simulations:int)-> None:   
    modules = [mod for mod in modules]
    os.makedirs(result_path, exist_ok=True)
    for mod in modules:
        mod_reuslt_path = os.path.join(result_path, mod.upper())
        os.makedirs(mod_reuslt_path, exist_ok=True)
        features_results = multi_time_simulation(tree_path,mod, features, seq_length=seq_length, indel_rate= indel_rate, module_res_path=mod_reuslt_path, number_of_simulations=number_of_simulations)
        save_results(features_results, mod_reuslt_path)