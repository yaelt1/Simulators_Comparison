import os
import subprocess

def AlisimCommandline(tree_path:str, result_path:str, seq_length:int=1000, indel_rate:float=0.01, indel_size = "POW{1.08/50},POW{1.08/50}",output_alignment_filename = 'alignment_w_indel'):
	"""
	Run AliSim
 
	:result_path = path to save results
 
	:seq_length = length of sequnces
	"""
	origin_dir = os.getcwd()
	os.chdir(result_path)
	print(os.getcwd())
	indel_rates =str(indel_rate)+","+str(indel_rate)
	run_iqtree_alisim(output_alignment_filename,tree_path, indel_rates, indel_size, seq_length=seq_length)
	


def run_iqtree_alisim(output_alignment_filename, tree, indel_rates, indel_size, seq_length:int=1000):
    command = f'iqtree2 --alisim {output_alignment_filename} -m JC -t {tree} --indel {indel_rates} --indel-size {indel_size} --out-format fasta --length {seq_length}'
    subprocess.run(command, shell=True)  