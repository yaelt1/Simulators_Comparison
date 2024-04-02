import os
import tempfile
import subprocess
from collections import OrderedDict
from indelible_dict import indelible_dict
import globals

def get_indelible_config():
	indelible_config = OrderedDict()
	indelible_config["[TYPE]"] = 'AMINOACID 2'
	indelible_config["[MODEL]"] = 'modelname'
	indelible_config["[submodel]"] = 'JTT'
	indelible_config["[indelmodel]"] = 'POW 1.08 50'
	indelible_config["[indelrate]"] = '0.01'
	indelible_config["[rates]"] = ' 0.25 0.50 10'
	indelible_config["[statefreq]"] = ' 0.25  0.25  0.25  0.25'
	indelible_config["[TREE]"] = 'treename (A:0.1,B:0.1);'
	indelible_config["[PARTITIONS]"] = 'partitionname\n[treename modelname 3000]'
	indelible_config["[EVOLVE]"] = "partitionname 100 outputname" + " " # note: indlible requires space at last command.
	return indelible_config


def prepare_indelible_control_file(result_path, model_parameters):
    indelible_config = get_indelible_config()

    indelible_config["[TREE]"] = f'treename {model_parameters["tree"]}'
    indelible_config["[PARTITIONS]"] = f'partitionname\n[treename modelname {model_parameters["length"]}]'
    indelible_config["[EVOLVE]"] = f'partitionname 1 outputname1' + " " 

    if model_parameters["mode"] == "amino":
        indelible_config["[TYPE]"] = 'AMINOACID 2'

    if ("w_indels" in model_parameters.keys()):
        indelible_config["[MODEL]"] = 'mymodel2'
        indelible_config["[submodel]"] = 'JTT'
        indelible_config["[indelmodel]"] = model_parameters["indelmodel"]
        indelible_config["[indelrate]"] = model_parameters["indelrate"]
        indelible_config["[PARTITIONS]"] = f'partitionname\n[treename mymodel2 {model_parameters["length"]}]'

    del indelible_config['[statefreq]']
    del indelible_config['[rates]']

    if model_parameters["mode"] == "nuc":
        indelible_config["[TYPE]"] = 'NUCLEOTIDE 2'
        
        if model_parameters["submodel"] == "GTR":
            gtr_params = ' '.join([f"{model_parameters['rates'][ind]:.9f}" for ind in range(5)])
            frequencies = ' '.join([f"{model_parameters['freq'][ind]:.6f}" for ind in range(4)])
            rates = f"{model_parameters['inv_prop']} {model_parameters['gamma_shape']} {model_parameters['gamma_cats']}"

            indelible_config["[submodel]"] = f'{model_parameters["submodel"]} {gtr_params}'
            indelible_config["[rates]"] = rates
            indelible_config["[statefreq]"] = frequencies

        if model_parameters["submodel"] == "JC":
            indelible_config["[submodel]"] = f'{model_parameters["submodel"]}'
            
            del indelible_config['[statefreq]']
            del indelible_config['[rates]']

    control_file_path = os.path.join(result_path, 'control.txt')
    os.makedirs(result_path, exist_ok=True)
    with open(control_file_path, 'w') as fout:
        for key in indelible_config:
            to_write = f'{key} {indelible_config[key]}\n'
            fout.write(to_write)

def IndelibleCommandline(result_path, model_params = indelible_dict):
	"""
	runs indelible.
	Requires control.txt at res_path and indelible command
	"""
	origin_dir = os.getcwd()
	prepare_indelible_control_file(result_path, model_params)
	os.chdir(result_path)
	subprocess.run(globals.INDELIBLE_EXECUTABLE, stdout=subprocess.DEVNULL, shell=True)
	indelible_msa_list = parse_indelible_output(result_path)
	return indelible_msa_list

def parse_indelible_output(res_path):
	"""
	reads the output of indelible and parse it to list of msas
	"""
	with open(os.path.join(res_path,'outputname1.fas'),'r') as f:
		indelible_subs = f.read()

	indelible_msa_list = [s[s.index("\n"):].replace("\n","") for s in indelible_subs.split(">")[1:]]
	
	return indelible_msa_list



def indelible_wrapper(result_path,tree_path, indel_rate: float = 0.01, seq_length: int = 1000, indel_dict=None):
    model_params = indel_dict if indel_dict is not None else indelible_dict
    with open(tree_path, 'r') as f:
        tree = f.read()
        model_params["tree"] = tree 
    model_params["indelrate"] = indel_rate
    model_params["length"] = seq_length
    result = IndelibleCommandline(model_params=model_params, result_path=result_path)
    return result