import pathlib
from itertools import repeat
import math
from _Sailfish import Tree, Simulator, SimProtocol, Msa, DiscreteDistribution, modelFactory, alphabetCode, modelCode
import random
import os, io, sys, tempfile
import time




def Sailfish_sim(tree_path:str,inse_r: int =0.01, del_r:int =0.01, n:int = 10000, seq_length:int = 1000,):
    tree = Tree(tree_path, True)
    root_node = tree.root
    protocol = SimProtocol(tree)

    rates_i = list(repeat(inse_r, tree.num_nodes))
    rates_d = list(repeat(del_r, tree.num_nodes))

    area_under = sum(i**-1.08 for i in range(1,51))

    length_probability_dist = [(i**-1.08)/area_under for i in range(1,51)]



    dist = DiscreteDistribution(length_probability_dist)
    dist.set_seed(random.randint(1, 1000))
    all_dists = [dist for i in range(tree.num_nodes)]
    protocol.set_sequence_size(seq_length)
    protocol.set_insertion_rates(rates_i)
    protocol.set_deletion_rates(rates_d)
    protocol.set_insertion_length_distributions(all_dists)
    protocol.set_deletion_length_distributions(all_dists)
    protocol.set_seed(random.randint(1, 100000))

  
    sim = Simulator(protocol)
    mFac = modelFactory(tree)

    mFac.set_alphabet(alphabetCode.AMINOACID)
    mFac.set_replacement_model(modelCode.JONES)
    mFac.set_gamma_parameters(0.5, 4)

    substitution_list = []

    blockmap = sim.gen_indels()
    msa = Msa(blockmap, root_node)
    sim.init_substitution_sim(mFac)
    substitutions = sim.gen_substitutions(msa.length())
    substitution_list.append(substitutions)
    
    msa.fill_substitutions(substitutions)
    
    msa.write_msa("msa.fas")
   



if __name__ == "__main__":
    Sailfish_sim()
