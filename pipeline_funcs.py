import sys
import os
import subprocess
import copy
import time
import random
import shutil
import numpy
from numpy.lib import recfunctions as rfn
import xml.etree.ElementTree as ET
import multiprocessing
from Bio import Phylo
from Bio import AlignIO
from datetime import datetime


def find_result_files(args, hypothesis_file_list):
	if args.ensemble_parts is not None and args.ensemble_parts >= 1:
		# Do analysis for ensemble model directory structure
		pass
		result_files = {}
		for hypothesis_file in hypothesis_file_list:
			hypothesis = os.path.splitext(os.path.basename(hypothesis_file))[0].replace("_hypothesis", "")
			result_files[hypothesis] = {}
			for i in range(1, args.ensemble_coverage+1):
				result_files[hypothesis][i] = {}
				for j in range(1, args.ensemble_parts+1):
					result_files[hypothesis][i][j] = {}
					result_files[hypothesis][i][j]["weights"] = os.path.join(args.output,
										     "{}_rep{}_part{}".format(args.output, i, j),
										     "{}_mapped_feature_weights.txt".format(hypothesis))
					result_files[hypothesis][i][j]["predictions"] = os.path.join(args.output,
										     "{}_rep{}_part{}".format(args.output, i, j),
										     "{}_gene_predictions.txt".format(hypothesis))
									     
	else:
		# Do analysis for singular model directory structure
		pass
		result_files = {}
		for hypothesis_file in hypothesis_file_list:
			hypothesis = os.path.splitext(os.path.basename(hypothesis_file))[0].replace("_hypothesis", "")
			result_files[hypothesis] = {}
			result_files[hypothesis]["weights"] = os.path.join(args.output,
									       "{}_mapped_feature_weights.txt".format(hypothesis))
			result_files[hypothesis]["predictions"] = os.path.join(args.output,
									       "{}_gene_predictions.txt".format(hypothesis))
	return result_files


def parse_result_files(args, file_dict):
	weights = {}
	for hypothesis in file_dict.keys():
		temp_weights = {}
		last_gene = ""
		for i in range(1, args.ensemble_coverage+1):
			for j in range(1, args.ensemble_parts+1):
				with open(file_dict[hypothesis][i][j]["weights"], 'r') as file:
					for line in file:
						data = line.strip().split("\t")
						rowname = data[0].split("_")
						feature = "_".join(rowname[-2:])
						gene = data[0].replace("_{}".format(feature), "")
						if i == 1:
							if gene != last_gene:
								last_gene = gene
								temp_weights[gene] = {}
							temp_weights[gene][feature] = []
						temp_weights[gene][feature].append(float(data[1]))
		weights[hypothesis] = temp_weights
	return weights


def analyze_ensemble_weights(args, weights):
	score_tables = {}
	for hypothesis in weights.keys():
		counts = {}
		totals = {}
		scores = {}
		outfile = os.path.join(args.output, "{}_weight_analysis.txt".format(hypothesis))
		for gene in weights[hypothesis].keys():
			counts[gene] = {}
			for feature in weights[hypothesis][gene].keys():
				counts[gene][feature] = sum([1 for x in weights[hypothesis][gene][feature] if abs(x)>0])
			totals[gene] = sum(counts[gene].values())
			scores[gene] = totals[gene]/len(counts[gene])
		with open(outfile, 'w') as file:
			file.write("{}\t{}\t{}\n".format("Gene","Total Non-zero","Score"))
			for gene in weights[hypothesis].keys():
				file.write("{}\t{}\t{}\n".format(gene, totals[gene], scores[gene]))
		score_tables[hypothesis] = {gene:(totals[gene], scores[gene]) for gene in weights[hypothesis].keys()}
	return score_tables


def apply_group_weights(aln_list_file, gene_sums, group_weights_file, species_list):
	group_sums = {}
	aln_list = []
	gene_list = list(gene_sums.keys())
	with open(aln_list_file, 'r') as file:
		for line in file:
			group_name = ",".join([os.path.splitext(os.path.basename(aln_filename))[0] for aln_filename in line.strip().split(",") if os.path.splitext(os.path.basename(aln_filename))[0] in gene_list])
			#if group_name != '':
			aln_list.append(group_name)
	group_weight_list = []
	with open(group_weights_file, 'r') as file:
		file.readline()
		file.readline()
		for weight in file.readline().strip().split(","):
			group_weight_list.append(float(weight))
	# print("aln_list len: {}\n group_weight_list len: {}".format(len(aln_list), len(group_weight_list)))
	for (aln_group, group_weight) in zip(aln_list, group_weight_list):
		if aln_group == '':
			continue
		group_sums[aln_group] = {}
		for seq_id in species_list:
			# group_sums[aln_group][seq_id] = sum([gene_sums[x].get(seq_id, numpy.nan) for x in aln_group.split(",")]) * group_weight
			group_sums[aln_group][seq_id] = sum([gene_sums[x].get(seq_id, 0) for x in aln_group.split(",")])
	return group_sums


def read_fasta(filename):
	with open(filename, 'r') as file:
		sequences = {}
		sequence = ''
		key = ''
		for line in file:
			if line.startswith('>'):
				if key and sequence:  # not empty
					sequences[key] = sequence
					sequence = ''  # reset for next sequence
				key = line.strip()[1:]  # remove '>' and newline character
			else:
				sequence += line.strip()  # concatenate lines of sequence
		if key and sequence:  # for the last sequence in file
			sequences[key] = sequence
	return sequences


def parse_response_file(response_filename, species_list):
	responses = {seq_id: None for seq_id in species_list}
	seq_order = []
	missing_seqids = []
	with open(response_filename, 'r') as file:
		custom_responses = [tuple(line.strip().split("\t")) for line in file]
		for custom_response in custom_responses:
			if custom_response[0] not in responses.keys():
				responses[custom_response[0]] = None
				missing_seqids = missing_seqids + [custom_response[0]]
			if responses[custom_response[0]] is None:
				responses[custom_response[0]] = custom_response[1]
				seq_order.append(custom_response[0])
			else:
				raise Exception("Response value of sequence {} specified more than once".format(custom_response[0]))
	# for seq_id in responses.keys():
	# 	if responses[seq_id] is None:
	# 		responses[seq_id] = "0"
	sorted_responses = {seq_id: responses[seq_id] for seq_id in seq_order}
	sorted_responses.update({seq_id: responses[seq_id] for seq_id in responses.keys() if seq_id not in seq_order})
	return sorted_responses, missing_seqids


def extract_gene_sums(aln_list, aln_lib, model):
	gene_files = {}
	gene_sums = {}
	gene_signifcance_scores = {}
	aln_dir = os.path.dirname(aln_list)
	with open(aln_list, 'r') as file:
		for line in file:
			for aln_filename in line.strip().split(","):
				gene_files[os.path.splitext(os.path.basename(aln_filename.strip()))[0]] = os.path.join(aln_dir, aln_filename.strip())
	found_gene_list = list(gene_files.keys())
	for gene in model.keys():
		if gene == "Intercept":
			continue
		elif gene not in found_gene_list:
			raise Exception("Gene {} present in model, but not present in input files.".format(gene))
	for gene in model.keys():
		if gene == "Intercept":
			continue
		gene_sums[gene] = {}
		if gene not in aln_lib.keys():
			aln_lib[gene] = read_fasta(gene_files[gene])
		for seq_id in aln_lib[gene].keys():
			gene_sums[gene][seq_id] = sum([model[gene][pos].get(aln_lib[gene][seq_id][pos], 0) for pos in model[gene].keys()])
		# gene_signifcance_scores[gene] = sum([sum(model[gene][pos].values()) for pos in model[gene].keys()])
		gene_signifcance_scores[gene] = sum([sum([abs(val) for val in pos.values()]) for pos in model[gene].values()])
	return gene_sums, gene_signifcance_scores


def read_ESL_model(filename):
	model = {}
	last_gene = ""
	last_pos = -1
	with open(filename, 'r') as file:
		for line in file:
			data = line.strip().split("\t")
			if data[0] == "Intercept":
				model["Intercept"] = float(data[1])
				continue
			feature = data[0].split("_")
			gene = "_".join(feature[0:-2])
			pos = int(feature[-2])
			allele = feature[-1]
			weight = float(data[1])
			if gene != last_gene:
				model[gene] = {pos: {allele: weight}}
			elif pos != last_pos:
				model[gene].update({pos: {allele: weight}})
			else:
				model[gene][pos].update({allele: weight})
			last_gene = gene
			last_pos = pos
	return model


def apply_ESL_model(aln_list, aln_lib, model_file, hypothesis_file, groups_filename, output_filename, missing_seqs):
	model_dir = os.path.dirname(model_file)
#	model_basename = os.path.splitext(os.path.basename(hypothesis_file))[0].replace("_hypothesis", "").replace("_mapped_feature_weights", "")
#	model = read_ESL_model(os.path.join(model_dir, model_basename + "_hypothesis_out_feature_weights.txt"))
	model = read_ESL_model(model_file)
	gene_sums, gene_significance_scores = extract_gene_sums(aln_list, aln_lib, model)
	species_list = set()
	for gene in gene_sums.keys():
		species_list.update(list(gene_sums[gene].keys()))
	response, missing_seqids = parse_response_file(hypothesis_file, species_list)
	for gene in gene_sums.keys():
		for missing_seqid in missing_seqids:
			gene_sums[gene][missing_seqid] = 0.0
	weighted_group_sums = apply_group_weights(aln_list, gene_sums, groups_filename, species_list)
	group_list = list(weighted_group_sums.keys())
	with open(output_filename, 'w') as file:
		# file.write("SeqID\tResponse\tPrediction\tIntercept\t{}\n".format("\t".join([",".join(x) for x in gene_list])))
		# file.write("SeqID\tResponse\tPrediction\tIntercept\t{}\n".format("\t".join(gene_list)))
		file.write("SeqID\tResponse\tPrediction\tIntercept\t{}\n".format("\t".join(group_list)))
		for seq_id in list(species_list) + missing_seqids:
			if response[seq_id] is None:
				continue
			prediction = model["Intercept"]
			for group in group_list:
				prediction += weighted_group_sums[group].get(seq_id, 0)
			#for i in range(0, len(gene_sums)):
			for group in group_list:
				if sum([1 for x in group.split(",") if (seq_id, x) in missing_seqs]) == len(group):
					#gene_sums[i] = numpy.nan
					weighted_group_sums[group] = numpy.nan
			file.write("{}\t{}\t{}\t{}\t{}\n".format(seq_id, response[seq_id], prediction, model["Intercept"],
													 "\t".join([str(weighted_group_sums[group].get(seq_id, numpy.nan)) for group in group_list])))
	if "gene_predictions_xval" not in output_filename:
		with open(str(output_filename).replace("_gene_predictions", "_GSS"), 'w') as file:
			file.write("{}\t{}\n".format("Gene","GSS"))
			# for (gene, weights) in zip(gene_list, group_weights):
			for gene in gene_significance_scores.keys():
				file.write("{}\t{}\n".format(gene, str(gene_significance_scores[gene])))


def generate_gene_prediction_table(weights_filename, responses_filename, groups_filename, features_filename, output_filename, gene_list, missing_seqs, group_list, model, features, field_filename=None):
	# Read weights, responses, and group indices files
	start = datetime.now()
	#model = xml_model_to_dict(weights_filename)
	weight_list = numpy.asarray(model["weight_list"])
	seqlist = []
	responses = {}
	with open(responses_filename, 'r') as file:
		for line in file:
			data = line.strip().split("\t")
			seqlist.append(data[0])
			responses[data[0]] = data[1]
	with open(groups_filename, 'r') as file:
		line1_data = [int(x) for x in file.readline().strip().split(",")]
		line2_data = [int(x) for x in file.readline().strip().split(",")]
		line3_data = [float(x) for x in file.readline().strip().split(",")]
	field = list(range(0, len(model["weight_list"])))
	if field_filename is not None:
		with open(field_filename, 'r') as file:
			field = [int(x)-1 for x in file.readline().strip().split(",") if len(x) > 0]
	group_indices = list(zip(line1_data, line2_data, line3_data))
	group_weights = []
	#print(len(model["weight_list"]))
	#print(len(field))
	for group in group_indices:
		#print(group)
		group_weights.append(numpy.asarray([model["weight_list"][field[x]] for x in range(group[0]-1, group[1])]))
	#features = numpy.loadtxt(features_filename, delimiter=',')
	# for (index, weight) in zip(group_indices, group_weights):
	# 	group_sums.append(numpy.sum(features[:, field[index[0] - 1:index[1]]] * weight, 1))
	#print(group_indices)
	group_sums = numpy.stack([numpy.sum(features[:, field[index[0] - 1:index[1]]] * weight, 1) for (index, weight) in zip(group_indices, group_weights)]).transpose()
	#predictions = numpy.sum(features * weight_list, 1) + model["intercept"]
	predictions = numpy.sum(group_sums, 1) + model["intercept"]
	# print("Time elapsed while multiplying features by group weights: {}".format(datetime.now() - start))
	with open(output_filename, 'w') as file:
		file.write("SeqID\tResponse\tPrediction\tIntercept\t{}\n".format("\t".join([",".join(x) for x in group_list])))
		for (seqid, gene_sums, prediction) in zip(seqlist, group_sums, predictions):
			for i in range(0, len(gene_sums)):
				if sum([1 for x in group_list[i] if (seqid, x) in missing_seqs]) == len(group_list[i]):
					#gene_sums[i] = "N/A"
					gene_sums[i] = numpy.nan
			file.write("{}\t{}\t{}\t{}\t{}\n".format(seqid, responses[seqid], prediction, model["intercept"], "\t".join([str(x) for x in gene_sums])))
	if "gene_predictions_xval" not in output_filename:
		with open(str(output_filename).replace("_gene_predictions", "_GSS"), 'w') as file:
			file.write("{}\t{}\n".format("Gene","GSS"))
			for (gene, weights) in zip(gene_list, group_weights):
				file.write("{}\t{}\n".format(gene, str(sum(numpy.abs(weights)))))


# Takes a list of alignment files and splits it into randomly selected subsets of equal size and returns the filenames
def split_gene_list(aln_list_filename, partitions):
	partition_files = []
	dirname = os.path.dirname(aln_list_filename)
	basename = "{}_{}".format(os.path.splitext(os.path.basename(aln_list_filename))[0], random.randrange(100000,999999))	
	with open(aln_list_filename, 'r') as file:
		aln_file_list = [] #List
		partitioned_file_list = [] #List of lists
		for line in file:
			aln_file_list.append(line.strip())
	random.shuffle(aln_file_list)
	partition_sizes = [0 for i in range(0, partitions)]
	for i in range(0, len(aln_file_list)):
		partition_sizes[i % partitions] += 1;
	for i in range(0, partitions):
		filename = "{}_partition_{}.txt".format(basename, i+1)
		filename = os.path.join(dirname, filename)
		with open(filename, 'w') as file:
			for j in range(0, partition_sizes[i]):
				file.write("{}\n".format(aln_file_list.pop()))
		partition_files.append(filename)
	return partition_files


def xml_model_to_dict(model_filename):
	params = {}
	params["intercept"] = 0
	params["lambda1"] = 0
	params["weight_list"] = []
	# Read weights and responses file
	xml_weights = ET.parse(model_filename)
	xml_tree = xml_weights.getroot()
	if xml_tree.tag != "model":
		raise Exception("Unexpected model XML format")
	for child1 in xml_tree:
		if child1.tag == "intercept_value":
			params["intercept"] = float(child1.text)
		elif child1.tag == "lambda1":
			params["lambda1"] = float(child1.text)
		elif child1.tag == "parameters":
			for child2 in child1:
				if child2.tag == "item":
					params["weight_list"].append(float(child2.text))
				else:
					params[child2.tag] = child2.text
	return params


def lookup_by_names(tree):
	names = {}
	for clade in tree.find_clades():
		if clade.name:
			if clade.name in names:
				raise ValueError("Duplicate key: %s" % clade.name)
			names[clade.name] = clade
	return names


def generate_hypothesis_set(args):
	newick_filename = args.tree
	nodelist_filename = args.nodelist
	response_filename = args.response
	auto_name_nodes = args.auto_name_nodes
	cladesize_cutoff_lower = args.cladesize_cutoff_lower
	cladesize_cutoff_upper = args.cladesize_cutoff_upper
	auto_name_length = args.auto_name_length
	smart_sampling = args.smart_sampling
	if args.slep_sample_balance or args.smart_sampling:
		slep_sample_balance = True
	else:
		slep_sample_balance = False
	tree = Phylo.parse(newick_filename, 'newick').__next__()
	taxa_list = [x.name for x in tree.get_terminals()]
	if cladesize_cutoff_upper is None:
		cladesize_cutoff_upper = len(taxa_list)
	taxa_list.reverse()
	auto_names = []
	if auto_name_nodes:
		i = 0
		for clade in tree.find_clades():
			if not clade.name:
				new_name = "{}_{}_{}".format(clade[0].get_terminals()[0].name[0:auto_name_length], clade[1].get_terminals()[0].name[0:auto_name_length], i)
				if new_name in auto_names:
					raise ValueError("Duplicate auto generated name: {}\nIncrease size of auto_name_length parameter and try again.".format(new_name))
				else:
					clade.name = new_name
					auto_names += [new_name]
					i += 1
		Phylo.write(tree, "auto_named_{}".format(os.path.basename(newick_filename)), "newick")
	nodes = lookup_by_names(tree)
	# print(tree)
	if nodelist_filename is None:
		nodelist = [key for key in nodes if key not in taxa_list]
	else:
		with open(nodelist_filename, 'r') as file:
			nodelist = [line.strip() for line in file]
	# print(nodelist)
	nodelist = [x for x in nodelist if
				len(nodes[x].get_terminals()) >= cladesize_cutoff_lower and len(nodes[x].get_terminals()) <= cladesize_cutoff_upper and len(
					nodes[x].get_terminals()) < len(taxa_list)]
	# print(nodelist)
	# print(tree)
	responses = {}
	if response_filename is None:
		distance_matrix = {t1: {t2: tree.distance(t1, t2) for t2 in taxa_list} for t1 in taxa_list}
		for nodename in nodelist:
			if smart_sampling is None:
				responses[nodename] = {x: -1 for x in taxa_list}
				for terminal in nodes[nodename].get_terminals():
					responses[nodename][terminal.name] = 1
			elif smart_sampling == 1:
				# print(nodename)
				responses[nodename] = {x: 0 for x in taxa_list}
				for terminal in nodes[nodename].get_terminals():
					responses[nodename][terminal.name] = 1
				target = nodes[nodename]
				response_sum = sum(responses[nodename].values())
				try:
					parent = tree.get_path(target)[-2]
					for cousin in parent.get_terminals():
						if responses[nodename][cousin.name] == 0:
							responses[nodename][cousin.name] = -1
				except:
					parent = tree.root
					for cousin in parent.get_terminals():
						if responses[nodename][cousin.name] == 0:
							responses[nodename][cousin.name] = -1
				response_sum = sum(responses[nodename].values())
				if response_sum < 0.1 * len(nodes[nodename].get_terminals()):
					temp_distance = copy.deepcopy(distance_matrix)
					negative_set = [key for key in responses[nodename].keys() if responses[nodename][key] == -1]
					for i in range(response_sum, 0):
						#get minimum pair distance from matrix
						#print(negative_set)
						#print(temp_distance.keys())
						minimum_distance = min([min([temp_distance[t1][t2] for t1 in negative_set if t1 != t2]) for t2 in negative_set])
						min_taxa = set()
						[[min_taxa.add(t1) for t1 in negative_set if temp_distance[t1][t2] == minimum_distance] for t2 in negative_set]
						min_taxa = list(min_taxa)
						#print(min_taxa)
						random.shuffle(min_taxa)
						#randomly delete half of pair
						#print(min_taxa[0])
						responses[nodename][min_taxa[0]] = 0
						negative_set.remove(min_taxa[0])
			elif smart_sampling == 2:
				# print(nodename)
				responses[nodename] = {x: 0 for x in taxa_list}
				for terminal in nodes[nodename].get_terminals():
					responses[nodename][terminal.name] = 1
				target = nodes[nodename]
				response_sum = sum(responses[nodename].values())
				while response_sum > 0.1 * len(nodes[nodename].get_terminals()):
					try:
						parent = tree.get_path(target)[-2]
						for cousin in parent.get_terminals():
							if responses[nodename][cousin.name] == 0:
								responses[nodename][cousin.name] = -1
						response_sum = sum(responses[nodename].values())
						target = parent
					except:
						parent = tree.root
						for cousin in parent.get_terminals():
							if responses[nodename][cousin.name] == 0:
								responses[nodename][cousin.name] = -1
						response_sum = 0
					if len(taxa_list) < 2.0 * len(nodes[nodename].get_terminals()):
						pass
				response_sum = sum(responses[nodename].values())
#				if response_sum < 0.1 * len(nodes[nodename].get_terminals()):
				if response_sum < 0:
					temp_distance = copy.deepcopy(distance_matrix)
					negative_set = [key for key in responses[nodename].keys() if responses[nodename][key] == -1]
					for i in range(response_sum, 0):
						#get minimum pair distance from matrix
						#print(negative_set)
						#print(temp_distance.keys())
						minimum_distance = min([min([temp_distance[t1][t2] for t1 in negative_set if t1 != t2]) for t2 in negative_set])
						min_taxa = set()
						[[min_taxa.add(t1) for t1 in negative_set if temp_distance[t1][t2] == minimum_distance] for t2 in negative_set]
						min_taxa = list(min_taxa)
						#print(min_taxa)
						random.shuffle(min_taxa)
						#randomly delete half of pair
						#print(min_taxa[0])
						responses[nodename][min_taxa[0]] = 0
						negative_set.remove(min_taxa[0])
				response_sum = sum(responses[nodename].values())
#				elif response_sum > 0.1 * len(nodes[nodename].get_terminals()):
				if response_sum > 0:
					temp_distance = copy.deepcopy(distance_matrix)
#					negative_set = [key for key in responses[nodename].keys() if responses[nodename][key] == -1]
					positive_set = [key for key in responses[nodename].keys() if responses[nodename][key] == 1]
					for i in range(0, response_sum):
						#get minimum pair distance from matrix
						#print(negative_set)
						#print(temp_distance.keys())
						minimum_distance = min([min([temp_distance[t1][t2] for t1 in positive_set if t1 != t2]) for t2 in positive_set])
						min_taxa = set()
						[[min_taxa.add(t1) for t1 in positive_set if temp_distance[t1][t2] == minimum_distance] for t2 in positive_set]
						min_taxa = list(min_taxa)
						#print(min_taxa)
						random.shuffle(min_taxa)
						#randomly delete half of pair
						#print(min_taxa[0])
						responses[nodename][min_taxa[0]] = 0
						positive_set.remove(min_taxa[0])
			elif smart_sampling == 3:
				# print(nodename)
				responses[nodename] = {x: 0 for x in taxa_list}
				for terminal in nodes[nodename].get_terminals():
					responses[nodename][terminal.name] = 1
				target = nodes[nodename]
				response_sum = sum(responses[nodename].values())
				while response_sum > 0.1 * len(nodes[nodename].get_terminals()):
					try:
						parent = tree.get_path(target)[-2]
						for cousin in parent.get_terminals():
							if responses[nodename][cousin.name] == 0:
								responses[nodename][cousin.name] = -1
						response_sum = sum(responses[nodename].values())
						target = parent
					except:
						parent = tree.root
						for cousin in parent.get_terminals():
							if responses[nodename][cousin.name] == 0:
								responses[nodename][cousin.name] = -1
						response_sum = 0
					if len(taxa_list) < 2.0 * len(nodes[nodename].get_terminals()):
						pass
				response_sum = sum(responses[nodename].values())
				if response_sum < 0.1 * len(nodes[nodename].get_terminals()):
					temp_distance = copy.deepcopy(distance_matrix)
					negative_set = [key for key in responses[nodename].keys() if responses[nodename][key] == -1]
					for i in range(response_sum, 0):
						#get minimum pair distance from matrix
						#print(negative_set)
						#print(temp_distance.keys())
						minimum_distance = min([min([temp_distance[t1][t2] for t1 in negative_set if t1 != t2]) for t2 in negative_set])
						min_taxa = set()
						[[min_taxa.add(t1) for t1 in negative_set if temp_distance[t1][t2] == minimum_distance] for t2 in negative_set]
						min_taxa = list(min_taxa)
						#print(min_taxa)
						random.shuffle(min_taxa)
						#randomly delete half of pair
						#print(min_taxa[0])
						responses[nodename][min_taxa[0]] = 0
						negative_set.remove(min_taxa[0])
	else:
		with open(response_filename, 'r') as file:
			basename = os.path.splitext(os.path.basename(response_filename))[0]
			responses[basename] = {x: None for x in taxa_list}
			custom_responses = [tuple(line.strip().split("\t")) for line in file]
			for response in custom_responses:
				for terminal in nodes[response[0]].get_terminals():
					if responses[basename][terminal.name] is None:
						responses[basename][terminal.name] = response[1]
					else:
						raise Exception("Response value of sequence {} specified more than once".format(terminal.name))
			for key in responses[basename].keys():
				if responses[basename][key] is None:
					responses[basename][key] = "0"
	hypothesis_file_list = []
	slep_opts_file_list = []
	for nodename in responses.keys():
		pos_idx = 1
		neg_idx = 1
		pos_idxs = []
		neg_idxs = []
		with open("{}_{}_hypothesis.txt".format(nodename, args.output), 'w') as file:
			for taxa in taxa_list:
				if responses[nodename][taxa] not in [0, "0"]:
					file.write("{}\t{}\n".format(taxa, responses[nodename][taxa]))
					if float(responses[nodename][taxa]) > 0:
						pos_idxs.append(pos_idx)
						pos_idx += 1
						if pos_idx > args.xval:
							pos_idx = 1
					elif float(responses[nodename][taxa]) < 0:
						neg_idxs.append(neg_idx)
						neg_idx += 1
						if neg_idx > args.xval:
							neg_idx = 1
		random.shuffle(pos_idxs)
		random.shuffle(neg_idxs)
		with open("{}_{}_xval_groups.txt".format(nodename, args.output), 'w') as file:
			for taxa in taxa_list:
				if responses[nodename][taxa] not in [0, "0"]:
					if float(responses[nodename][taxa]) > 0:
						#file.write("{}\t{}\n".format(taxa, pos_idxs.pop()))
						file.write("{}\n".format(pos_idxs.pop()))
					elif float(responses[nodename][taxa]) < 0:
						#file.write("{}\t{}\n".format(taxa, neg_idxs.pop()))
						file.write("{}\n".format(neg_idxs.pop()))
		hypothesis_file_list += ["{}_{}_hypothesis.txt".format(nodename, args.output)]
		if slep_sample_balance:
			slep_opts_file_list += ["{}_{}_slep_opts.txt".format(nodename, args.output)]
			with open("{}_{}_slep_opts.txt".format(nodename, args.output), 'w') as opts_file:
				if args.slep_opts is not None:
					with open(args.slep_opts, 'r') as base_opts_file:
						for line in base_opts_file:
							opts_file.write(line)
				# ratio = sum([1.0 for x in responses[nodename].values() if x == 1])/sum([1.0 for x in responses[nodename].values() if x == -1])
				opts_file.write("{}_{}_sweights.txt\n".format(nodename, args.output))
				with open("{}_{}_sweights.txt".format(nodename, args.output), 'w') as sweights_file:
					sweights_file.write("{}\n".format(sum([1.0 for x in responses[nodename].values() if int(x) == 1])/sum([1.0 for x in responses[nodename].values() if int(x) == -1])))
					sweights_file.write("{}\n".format(1.0))
		else:
			slep_opts_file_list += [args.slep_opts]
	return hypothesis_file_list, slep_opts_file_list


def split_path(path):
	path_list = []
	while 1:
		path, dir = os.path.split(path)
		path_list.append(dir)
		if path == "":
			break
	path_list.reverse()
	return path_list


def generate_input_matrices(alnlist_filename, hypothesis_filename_list, args):
	output_basename = args.output
	options = ""
	modified_response = True
	if args.upsample_balance:
		options = "{} {}".format(options.strip(),"ub")
		modified_response = True
	elif args.downsample_balance:
		options = "{} {}".format(options.strip(),"db")
		modified_response = True
	if args.fuzz_indels:
		options = "{} {}".format(options.strip(),"fuzzIndels")
	response_file_list = []
	group_indices_file_list = []
	features_file_list = []
	aln_file_list = {}
	gene_list = []
	group_list = []
	field_file_list = []
	# Generate gene list from alignment list file
	with open(alnlist_filename) as file:
		for line in file:
			group = []
			for aln_filename in line.strip().split(","):
				basename = os.path.splitext(os.path.basename(aln_filename.strip()))[0]
				group.append(basename)
				if basename in aln_file_list.keys():
					if aln_file_list[basename] != aln_filename.strip():
						raise Exception("Found multiple alignment files with identical basename {}.".format(basename))
				else:
					aln_file_list[basename] = aln_filename.strip()
			gene_list.extend(group)
			group_list.append(group)
	preprocess_exe = os.path.join(os.getcwd(), "bin", "preprocess")
	preprocess_cwd, alnlist_filename = os.path.split(alnlist_filename)
	if preprocess_cwd == "":
			preprocess_cwd = "."
	if not modified_response:
		# Construct preprocessing command for first hypothesis file
		preprocess_cmd = "{}*{}*{}*{}*{}".format(preprocess_exe, os.path.join(os.getcwd(), hypothesis_filename_list[0]), alnlist_filename, output_basename, options)
		print(preprocess_cmd.replace("*"," "))
		subprocess.call(preprocess_cmd.split("*"), stderr=subprocess.STDOUT, cwd=preprocess_cwd)
		# Move generated inputs to top level directory
		if preprocess_cwd != ".":
			shutil.move(os.path.join(preprocess_cwd, output_basename), ".")
		# Substitute gene/group penalties
		if args.gene_penalties is not None:
			gene_penalties = {}
			penalty_list = []
			with open(os.path.join(output_basename, "group_indices_" + output_basename + ".txt"), 'r') as groups_file:
				lines = [val.strip() for val in groups_file.readlines()]
			with open(args.gene_penalties, 'r') as penalties_file:
				for line in penalties_file.readlines():
					data = line.strip().split('\t')
					gene_penalties[data[0]] = float(data[1])
			with open(os.path.join(preprocess_cwd, alnlist_filename), 'r') as file:
				for line in file.readlines():
					try:
						penalty_list.append(sum([gene_penalties[os.path.splitext(os.path.basename(filename))[0]] for filename in line.strip().split(',')]))
					except:
						raise Exception("Could not find gene penalty for one of the following genes: {}".format(','.join([os.path.splitext(os.path.basename(filename))[0] for filename in line.strip().split(',')])))
			with open(os.path.join(output_basename, "group_indices_" + output_basename + ".txt"), 'w') as groups_file:
				groups_file.write("{}\n".format(lines[0]))
				groups_file.write("{}\n".format(lines[1]))
				groups_file.write("{}\n".format('\t'.join([str(x) for x in penalty_list])))
		# Construct response input file for each additional hypothesis file
		for filename in hypothesis_filename_list:
			with open(filename, 'r') as infile:
				temp_fname = os.path.join(output_basename, "response_" + os.path.splitext(os.path.basename(filename))[0] + ".txt")
				with open(temp_fname, 'w') as outfile:
					for line in infile:
						outfile.write("{}\n".format(line.strip().split("\t")[1]))
					response_file_list.append(temp_fname)
					group_indices_file_list.append(os.path.join(output_basename, "group_indices_" + output_basename + ".txt"))
					features_file_list.append(os.path.join(output_basename, "feature_" + output_basename + ".txt"))
					field_file_list.append(os.path.join(output_basename, "field_" + output_basename + ".txt"))
		#return [os.path.join(output_basename, "feature_" + output_basename + ".txt"), os.path.join(output_basename, "group_indices_" + output_basename + ".txt"), response_file_list, gene_list]
		return [features_file_list, group_indices_file_list, response_file_list, gene_list, field_file_list, group_list]
	else:
		for filename in hypothesis_filename_list:
			# Construct preprocessing command
			preprocess_cmd = "{}*{}*{}*{}*{}".format(preprocess_exe, os.path.join(os.getcwd(), filename), alnlist_filename, output_basename, options)
			print(preprocess_cmd.replace("*"," "))
			hypothesis_basename = os.path.splitext(os.path.basename(filename))[0]
			if args.skip_preprocessing:
				if os.path.exists(os.path.join(output_basename, "feature_" + hypothesis_basename + ".txt")):
					print("Features file detected, skipping preprocessing step...")
				else:
					raise Exception("Preprocessing skipped, but no features file detected at {}.".format(os.path.join(output_basename, "feature_" + hypothesis_basename + ".txt")))
			else:
				subprocess.call(preprocess_cmd.split("*"), stderr=subprocess.STDOUT, cwd=preprocess_cwd)
				if preprocess_cwd != ".":
					if os.path.exists(output_basename):
						shutil.copytree(os.path.join(preprocess_cwd, output_basename), output_basename, dirs_exist_ok=True)
						shutil.rmtree(os.path.join(preprocess_cwd, output_basename))
					else:
						shutil.move(os.path.join(preprocess_cwd, output_basename), output_basename)
				if True:
					position_stats = {}
					#stat_keys = ["mic", "entropy"]
					stat_keys = ["mic"]
					for aln_basename in aln_file_list.keys():
						position_stats[aln_basename] = calculate_position_stats(os.path.join(preprocess_cwd, aln_file_list[aln_basename]), filename)
					with open(os.path.join(output_basename, "pos_stats_" + hypothesis_basename + ".txt"), 'w') as file:
						file.write("{}\t{}\n".format("Position Name", '\t'.join(stat_keys)))
						for aln_basename in position_stats.keys():
							for i in range(0, len(position_stats[aln_basename]["mic"])):
								file.write("{}_{}\t{}\n".format(aln_basename, i, '\t'.join([str(position_stats[aln_basename][stat_key][i]) for stat_key in stat_keys])))
				shutil.move(os.path.join(output_basename, "feature_" + output_basename + ".txt"), os.path.join(output_basename, "feature_" + hypothesis_basename + ".txt"))
				shutil.move(os.path.join(output_basename, "group_indices_" + output_basename + ".txt"), os.path.join(output_basename, "group_indices_" + hypothesis_basename + ".txt"))
				shutil.move(os.path.join(output_basename, "response_" + output_basename + ".txt"), os.path.join(output_basename, "response_" + hypothesis_basename + ".txt"))
				shutil.move(hypothesis_basename.replace("hypothesis", "xval_groups.txt"), os.path.join(output_basename, hypothesis_basename.replace("hypothesis", "xval_groups.txt")))
				shutil.move(os.path.join(output_basename, "field_" + output_basename + ".txt"), os.path.join(output_basename, "field_" + hypothesis_basename + ".txt"))
				shutil.move(os.path.join(output_basename, "feature_mapping_" + output_basename + ".txt"), os.path.join(output_basename, "feature_mapping_" + hypothesis_basename + ".txt"))
				try:
					shutil.move(os.path.join(output_basename, "resampled_" + output_basename + ".txt"), os.path.join(output_basename, "resampled_" + hypothesis_basename + ".txt"))
				except:
					pass
			response_file_list.append(os.path.join(output_basename, "response_" + hypothesis_basename + ".txt"))
			group_indices_file_list.append(os.path.join(output_basename, "group_indices_" + hypothesis_basename + ".txt"))
			field_file_list.append(os.path.join(output_basename, "field_" + hypothesis_basename + ".txt"))
			features_file_list.append(os.path.join(output_basename, "feature_" + hypothesis_basename + ".txt"))
		return [features_file_list, group_indices_file_list, response_file_list, gene_list, field_file_list, group_list]


def calculate_position_stats(aln_filename, hypothesis_filename):
	mic = []
	entropy = []
	with open(hypothesis_filename) as hypothesis_file:
		responses = {}
		species_counts = [0, 0]
		for line in hypothesis_file.readlines():
			data = line.strip().split('\t')
			if int(data[1]) in [-1, 1]:
				responses[data[0]] = int(data[1])
				if int(data[1]) == 1:
					species_counts[0] += 1
				elif int(data[1]) == -1:
					species_counts[1] += 1
			else:
				responses[data[0]] = 0
		aln = AlignIO.read(aln_filename, "fasta")
		aln.sort(key=lambda record: responses.get(record.id, 0))
		aln1 = aln[0:species_counts[0]]
		aln.sort(key=lambda record: responses.get(record.id, 0), reverse=True)
		aln2 = aln[0:species_counts[1]]
		for i in range(0, aln.get_alignment_length()):
			base_counts1 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
			base_counts2 = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
			for base in aln1[:, i]:
				base_counts1[base] = base_counts1.get(base, 0) + 1
			for base in aln2[:, i]:
				base_counts2[base] = base_counts2.get(base, 0) + 1
			base_counts1['total'] = sum([base_counts1[base] for base in ["A", "T", "C", "G"]])
			base_counts2['total'] = sum([base_counts2[base] for base in ["A", "T", "C", "G"]])
			base_counts = {base: base_counts1[base]+base_counts2[base] for base in ['A', 'T', 'C', 'G', 'total']}
			for counts in [base_counts1, base_counts2, base_counts]:
				for base in ['A', 'T', 'C', 'G']:
					if counts['total'] > 0:
						counts[base] = counts[base]/counts['total']
					else:
						counts[base] = 0
			ic1 = sum([0 if base_counts1[base] == 0 else base_counts1[base] * numpy.log10(base_counts1[base]) for base in ['A', 'T', 'C', 'G']])
			ic2 = sum([0 if base_counts2[base] == 0 else base_counts2[base] * numpy.log10(base_counts2[base]) for base in ['A', 'T', 'C', 'G']])
			ic = sum([0 if base_counts[base] == 0 else base_counts[base] * numpy.log10(base_counts[base]) for base in ['A', 'T', 'C', 'G']])
			mic.append(max(ic - ic1 - ic2, 0))
	return {"mic": mic, "entropy": entropy}


def run_esl(features_filename_list, groups_filename_list, response_filename_list, field_filename_list, args, slep_opts_filename_list, z_ind = None, y_ind = None):
	sparsity = args.lambda1
	group_sparsity = args.lambda2
	method = args.method
	if method == "leastr":
		method = "sg_lasso_leastr"
	elif method == "logistic":
		method = "sg_lasso"
	elif method == "ol_leastr":
		method = "overlapping_sg_lasso_leastr"
	elif method == "ol_logistic":
		method = "overlapping_sg_lasso_logisticr"
	else:
		raise Exception("Provided method name not recognized, please provide a valid method name.")
	weights_file_list = []
	esl_exe = os.path.join(os.getcwd(), "bin", method)
	# Run sg_lasso for each response file in response_filename_list
	for response_filename, features_filename, groups_filename, field_filename, slep_opts_filename in zip(response_filename_list, features_filename_list, groups_filename_list, field_filename_list, slep_opts_filename_list):
		if z_ind is None or y_ind is None:
			basename = str(os.path.splitext(os.path.basename(response_filename))[0]).replace("response_","")
		else:
			basename = str(os.path.splitext(os.path.basename(response_filename))[0]).replace("response_", "") + "_{}_{}".format(z_ind, y_ind)
		if slep_opts_filename is None:
			esl_cmd = "{}*-f*{}*-z*{}*-y*{}*-n*{}*-r*{}*-w*{}".format(esl_exe, features_filename, sparsity, group_sparsity, groups_filename, response_filename, basename + "_out_feature_weights")
		else:
			esl_cmd = "{}*-f*{}*-z*{}*-y*{}*-n*{}*-r*{}*-s*{}*-w*{}".format(esl_exe, features_filename, sparsity, group_sparsity, groups_filename, response_filename, slep_opts_filename, basename + "_out_feature_weights")
		if method in ["overlapping_sg_lasso_leastr", "overlapping_sg_lasso_logisticr"]:
			esl_cmd = esl_cmd + "*-g*{}".format(field_filename)
		if args.xval > 1:
			esl_cmd = esl_cmd + "*-x*{}".format(response_filename.replace("response_", "").replace("hypothesis.txt", "xval_groups.txt"))
		esl_cmd = esl_cmd + "*--model_format*flat"
		print(esl_cmd.replace("*"," "))
#		subprocess.call("touch {}".format(basename + "_out_feature_weights.xml"), stderr=subprocess.STDOUT, shell=True)
		subprocess.call(esl_cmd.split("*"), stderr=subprocess.STDOUT)
		weights_file_list.append(basename + "_out_feature_weights.txt")
	return weights_file_list


def run_esl_grid(features_filename_list, groups_filename_list, response_filename_list, field_filename_list, sparsity, group_sparsity, method, slep_opts_filename_list, lambda_list_filename, args, z_ind = None, y_ind = None):
	if method == "leastr":
		method = "sg_lasso_leastr"
	elif method == "logistic":
		method = "sg_lasso"
	elif method == "ol_leastr":
		method = "overlapping_sg_lasso_leastr"
	elif method == "ol_logistic":
		method = "overlapping_sg_lasso_logisticr"
	else:
		raise Exception("Provided method name not recognized, please provide a valid method name.")
	weights_file_list = []
	esl_exe = os.path.join(os.getcwd(), "bin", method)
	lambda_list = []
	with open(lambda_list_filename, 'r') as file:
		for line in file:
			data = line.strip().split("\t")
			lambda_list.append((data[0], data[1]))
	# Run sg_lasso for each response file in response_filename_list
	for response_filename, features_filename, groups_filename, field_filename, slep_opts_filename in zip(response_filename_list, features_filename_list, groups_filename_list, field_filename_list, slep_opts_filename_list):
		if z_ind is None or y_ind is None:
			basename = str(os.path.splitext(os.path.basename(response_filename))[0]).replace("response_","")
		else:
			basename = str(os.path.splitext(os.path.basename(response_filename))[0]).replace("response_", "") + "_{}_{}".format(z_ind, y_ind)
		if slep_opts_filename is None:
			esl_cmd = "{}*-f*{}*-z*{}*-y*{}*-n*{}*-r*{}*-w*{}".format(esl_exe, features_filename, sparsity, group_sparsity, groups_filename, response_filename, basename + "_out_feature_weights")
		else:
			esl_cmd = "{}*-f*{}*-z*{}*-y*{}*-n*{}*-r*{}*-s*{}*-w*{}".format(esl_exe, features_filename, sparsity, group_sparsity, groups_filename, response_filename, slep_opts_filename, basename + "_out_feature_weights")
		if method in ["overlapping_sg_lasso_leastr", "overlapping_sg_lasso_logisticr"]:
			esl_cmd = esl_cmd + "*-g*{}".format(field_filename)
		esl_cmd = esl_cmd + "*-l*{}".format(lambda_list_filename)
		if args.grid_gene_threshold:
			esl_cmd = esl_cmd + "*-c*{}".format(args.grid_gene_threshold)
		esl_cmd = esl_cmd + "*--model_format*flat"
		print(esl_cmd.replace("*"," "))
		#subprocess.call("touch {}".format(basename + "_out_feature_weights.xml"), stderr=subprocess.STDOUT, shell=True)
		if args.skip_processing:
			for model_file in ["{}_out_feature_weights_{}_{}.txt".format(basename, lambda_val2label(val[0]), lambda_val2label(val[1])) for val in lambda_list]:
				if not os.path.exists(model_file):
					raise Exception("Processing skipped, but no model file detected at {}.".format(model_file))
			print("All expected model files detected, skipping processing step...")
		else:
			subprocess.call(esl_cmd.split("*"), stderr=subprocess.STDOUT)
		weights_file_list.append(["{}_out_feature_weights_{}_{}.txt".format(basename, lambda_val2label(val[0]), lambda_val2label(val[1])) for val in lambda_list])
	return weights_file_list


def lambda_val2label(lambda_val):
	lambda_val = float(lambda_val)
	if "{:g}".format(lambda_val)[0:2] == "0.":
		return "{:g}".format(lambda_val)[2:]
	else:
		return "{:g}".format(lambda_val)


def process_grid_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs, group_list, args):
	missing_results = []
	aln_lib = {}
	for (weights_filename_list, hypothesis_filename, groups_filename, features_filename) in zip(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list):
		# outname = hypothesis_filename.replace("_hypothesis.txt", "_gene_predictions.txt")
		start = datetime.now()
		features = numpy.loadtxt(features_filename, delimiter=',')
		print("Time elapsed loading features file: {}".format(datetime.now() - start))
		for weights_filename in weights_filename_list:
			outname = weights_filename.replace("_hypothesis", "").replace("_out_feature_weights", "_gene_predictions").replace(".xml", ".txt")
			# generate_gene_prediction_table(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, groups_filename.replace("group_indices_", "field_"))
			# total_significance = generate_mapped_weights_file(weights_filename, groups_filename.replace("group_indices_", "feature_mapping_"))
			if os.path.exists(weights_filename):
				HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + process_single_grid_weight(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, features, args)
			elif args.grid_gene_threshold is not None:
				print("No results file detected, most likely due to grid_gene_threshold: {}".format(weights_filename))
				missing_results.append(weights_filename)
			else:
				raise Exception("Missing results file: {}".format(weights_filename))
	#	pool_results[hypothesis_filename].append(thread_pool.apply_async(process_single_grid_weight, args=(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, args)))
	# thread_pool.close()
	# thread_pool.join()
	# for hypothesis_filename in hypothesis_file_list:
	# 	for pool_result in pool_results[hypothesis_filename]:
	# 		HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + pool_result.get()
	return missing_results


def process_single_grid_weight(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, features, args):
	start = datetime.now()
	model = xml_model_to_dict(weights_filename)
	generate_gene_prediction_table(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, model, features, groups_filename.replace("group_indices_", "field_"))
	print("Time elapsed while generating gene prediction table: {}".format(datetime.now() - start))
	start = datetime.now()
	# total_significance = generate_mapped_weights_file(weights_filename, groups_filename.replace("group_indices_", "feature_mapping_"), model)
	total_significance = generate_significance_scores(weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.xml", "_mapped_feature_weights.txt"))
	print("Time elapsed while generating mapped weights file: {}".format(datetime.now() - start))
	if not args.preserve_xml:
		os.remove(weights_filename)
	return total_significance

# def           process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs, group_list):
def process_sparse_grid_weights(weights_file_list, hypothesis_file_list, groups_filename_list, HSS, missing_seqs, args):
	missing_results = []
	aln_lib = {}
	#for (weights_filename_list, hypothesis_filename, groups_filename, features_filename) in zip(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list):
	for (weights_filename_list, hypothesis_filename, groups_filename) in zip(weights_file_list, hypothesis_file_list, groups_filename_list):
		# outname = hypothesis_filename.replace("_hypothesis.txt", "_gene_predictions.txt")
		for weights_filename in weights_filename_list:
			outname = weights_filename.replace("_hypothesis", "").replace("_out_feature_weights", "_gene_predictions").replace(".xml", ".txt")
			# generate_gene_prediction_table(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, groups_filename.replace("group_indices_", "field_"))
			# total_significance = generate_mapped_weights_file(weights_filename, groups_filename.replace("group_indices_", "feature_mapping_"))
			if os.path.exists(weights_filename):
#				HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + process_single_grid_weight(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, features, args)
				HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + process_sparse_single_grid_weight(weights_filename, hypothesis_filename, groups_filename, outname, missing_seqs, aln_lib, args)
			elif args.grid_gene_threshold is not None:
				print("No results file detected, most likely due to grid_gene_threshold: {}".format(weights_filename))
				missing_results.append(weights_filename)
			else:
				raise Exception("Missing results file: {}".format(weights_filename))
	#	pool_results[hypothesis_filename].append(thread_pool.apply_async(process_single_grid_weight, args=(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, args)))
	# thread_pool.close()
	# thread_pool.join()
	# for hypothesis_filename in hypothesis_file_list:
	# 	for pool_result in pool_results[hypothesis_filename]:
	# 		HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + pool_result.get()
	return missing_results


def process_sparse_single_grid_weight(weights_filename, hypothesis_filename, groups_filename, outname, missing_seqs, aln_lib, args):
#	start = datetime.now()
#	model = xml_model_to_dict(weights_filename)
	apply_ESL_model(args.aln_list, aln_lib, weights_filename, hypothesis_filename, groups_filename, outname, missing_seqs)
#	generate_gene_prediction_table(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, model, features, groups_filename.replace("group_indices_", "field_"))
#	print("Time elapsed while generating gene prediction table: {}".format(datetime.now() - start))
#	start = datetime.now()
	# total_significance = generate_mapped_weights_file(weights_filename, groups_filename.replace("group_indices_", "feature_mapping_"), model)
	#total_significance = generate_significance_scores(weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.", "_mapped_feature_weights."))
	total_significance = generate_significance_scores(weights_filename, groups_filename)
#	print("Time elapsed while generating mapped weights file: {}".format(datetime.now() - start))
#	if not args.preserve_xml:
#		os.remove(weights_filename)
	return total_significance


def process_sparse_weights(weights_file_list, hypothesis_file_list, groups_filename_list, HSS, missing_seqs, args):
	aln_list = args.aln_list
	aln_lib = {}
	for (weights_filename, hypothesis_filename, groups_filename) in zip(weights_file_list, hypothesis_file_list, groups_filename_list):
		outname = weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.txt", "_gene_predictions.txt")
		apply_ESL_model(aln_list, aln_lib, weights_filename, hypothesis_filename, groups_filename, outname, missing_seqs)
		# total_significance = generate_significance_scores(weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.xml", "_mapped_feature_weights.txt"))
		total_significance = generate_significance_scores(weights_filename, groups_filename)
		shutil.move(weights_filename.replace("_hypothesis_out_feature_weights", "_PSS"), args.output)
		#os.remove(weights_filename)
		shutil.move(weights_filename, os.path.join(os.path.dirname(groups_filename), weights_filename.replace("_hypothesis_out_feature_weights.xml", "_mapped_feature_weights.txt")))
		HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + total_significance


def process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs, group_list):
	aln_lib = {}
	for (weights_filename, hypothesis_filename, groups_filename, features_filename) in zip(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list):
		# outname = hypothesis_filename.replace("_hypothesis.txt", "_gene_predictions.txt")
		outname = weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.xml", "_gene_predictions.txt")
		model = xml_model_to_dict(weights_filename)
		features = numpy.loadtxt(features_filename, delimiter=',')
		generate_gene_prediction_table(weights_filename, hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, model, features, groups_filename.replace("group_indices_", "field_"))
		# apply_ESL_model(aln_list, aln_lib, model_file, hypothesis_filename, outname, missing_seqs)
		total_significance = generate_mapped_weights_file(weights_filename, groups_filename.replace("group_indices_", "feature_mapping_"), model)
		# total_significance = generate_significance_scores(weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.xml", "_mapped_feature_weights.txt"))
		os.remove(weights_filename)
		HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + total_significance


def process_sparse_xval_weights(weights_file_list, hypothesis_file_list, groups_filename_list, aln_list,                 gene_list, xval, missing_seqs, group_list):
	aln_lib = {}
#	for (weights_filename, hypothesis_filename, groups_filename, features_filename) in zip(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list):
	for (weights_filename, hypothesis_filename, groups_filename) in zip(weights_file_list, hypothesis_file_list, groups_filename_list):
		# outname = hypothesis_filename.replace("_hypothesis.txt", "_gene_predictions.txt")
		# features = numpy.loadtxt(features_filename, delimiter=',')
		groups = numpy.loadtxt(hypothesis_filename.replace("_hypothesis", "_xval_groups"))
		xval_predictions = []
		header = ""
		for xval_id in range(1, xval+1):
			outname = weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.txt", "_gene_predictions_xval_{}.txt".format(xval_id))
			xval_hypothesis_filename = hypothesis_filename.replace("hypothesis.txt", "hypothesis_xval_{}.txt".format(xval_id))
			xval_model_filename = weights_filename.replace("_out_feature_weights.txt", "_out_feature_weights_xval_{}.txt".format(xval_id))
			with open(hypothesis_filename, 'r') as in_file:
				with open(xval_hypothesis_filename, 'w') as out_file:
					for i in range(0, len(groups)):
						line = in_file.readline()
						if groups[i] == xval_id:
							out_file.write(line)
			apply_ESL_model(aln_list, aln_lib, xval_model_filename, xval_hypothesis_filename, groups_filename, outname, missing_seqs)
			#total_significance = generate_mapped_weights_file(weights_filename, groups_filename.replace("group_indices_", "feature_mapping_"), model)
			os.remove(weights_filename.replace("_out_feature_weights.txt", "_out_feature_weights_xval_{}.txt".format(xval_id)))
			xval_predictions.append(numpy.genfromtxt(outname, delimiter='	', names=True, dtype=None, encoding=None))
			if xval_id == 1:
				with open(outname, 'r') as file:
					header = file.readline().strip()
			os.remove(outname)
			os.remove(xval_hypothesis_filename)
		#HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + total_significance
		xval_predictions_fname = weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.txt", "_gene_predictions_xval.txt")
		all_genes = set()
		for prediction in xval_predictions:
			all_genes.update(prediction.dtype.names)
		all_columns = ['SeqID', 'Response', 'Prediction', 'Intercept']
		for gene in all_genes:
			if gene not in ['SeqID', 'Response', 'Prediction', 'Intercept']:
				all_columns.append(gene)
		combined_prediction = numpy.concatenate([add_missing_columns(xval_prediction, all_columns) for xval_prediction in xval_predictions])
		with open(xval_predictions_fname, 'w') as file:
			file.write(header + '\n')
			numpy.savetxt(file, combined_prediction, delimiter='	', fmt='%s')


# Function to add missing columns with zeros
def add_missing_columns(data, columns):
	missing_columns = set(columns) - set(data.dtype.names)
	for column in missing_columns:
		data = rfn.append_fields(data, column, [0]*len(data), dtypes=None, usemask=False)
	data_reordered = numpy.empty(len(data), dtype=[(columns[0], 'U128')] + [(colname,'f8') for colname in columns if colname != 'SeqID'])
	for column_name in columns:
		data_reordered[column_name] = data[column_name]
	return data_reordered


def process_xval_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, xval, missing_seqs, group_list):
	for (weights_filename, hypothesis_filename, groups_filename, features_filename) in zip(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list):
		# outname = hypothesis_filename.replace("_hypothesis.txt", "_gene_predictions.txt")
		features = numpy.loadtxt(features_filename, delimiter=',')
		groups = numpy.loadtxt(hypothesis_filename.replace("_hypothesis", "_xval_groups"))
		xval_predictions = []
		header = ""
		for xval_id in range(1, xval+1):
			outname = weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.xml", "_gene_predictions_xval_{}.txt".format(xval_id))
			xval_hypothesis_filename = hypothesis_filename.replace("hypothesis.txt", "hypothesis_xval_{}.txt".format(xval_id))
			with open(hypothesis_filename, 'r') as in_file:
				with open(xval_hypothesis_filename, 'w') as out_file:
					for i in range(0, len(groups)):
						line = in_file.readline()
						if groups[i] == xval_id:
							out_file.write(line)
			model = xml_model_to_dict(weights_filename.replace("_out_feature_weights.xml", "_out_feature_weights_xval_{}.xml".format(xval_id)))
			generate_gene_prediction_table(weights_filename, xval_hypothesis_filename, groups_filename, features_filename, outname, gene_list, missing_seqs, group_list, model, features[numpy.where(groups == xval_id)[0]], groups_filename.replace("group_indices_", "field_"))
			#total_significance = generate_mapped_weights_file(weights_filename, groups_filename.replace("group_indices_", "feature_mapping_"), model)
			os.remove(weights_filename.replace("_out_feature_weights.xml", "_out_feature_weights_xval_{}.xml".format(xval_id)))
			xval_predictions.append(numpy.loadtxt(outname, delimiter='	', skiprows=1, dtype=str))
			if xval_id == 1:
				with open(outname, 'r') as file:
					header = file.readline().strip()
			os.remove(outname)
			os.remove(xval_hypothesis_filename)
		#HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + total_significance
		xval_predictions_fname = weights_filename.replace("_hypothesis", "").replace("_out_feature_weights.xml", "_gene_predictions_xval.txt")
		with open(xval_predictions_fname, 'w') as file:
			file.write(header + '\n')
			numpy.savetxt(file, numpy.vstack(xval_predictions), delimiter='	', fmt='%s')


def generate_significance_scores(model_filename, groups_filename):
	# Read weights and feature mapping files
	PSS = {}
	posname_list = []
	last_posname = ""
	feature_map = {}
	pos_stats = {}
	# output_filename = str(weights_filename).replace("_hypothesis_out_feature_weights.xml", "_mapped_feature_weights.txt")
#	output_filename = str(groups_filename).replace("_hypothesis.txt", "_mapped_feature_weights.txt").replace("group_indices_", "")
	with open(model_filename, 'r') as file:
		for line in file:
			data = line.strip().split("\t")
			if len(data) == 2:
				feature_map[data[0]] = float(data[1])
	if os.path.exists(groups_filename.replace("group_indices_", "pos_stats_")):
		with open(groups_filename.replace("group_indices_", "pos_stats_"), 'r') as file:
			for line in file:
				data = line.strip().split("\t")
				if len(data) > 1:
					pos_stats[data[0]] = data[1:]
	for feature in feature_map.keys():
		if feature == "Intercept":
			continue
		posname = feature[0:-2]
		if posname != last_posname:
			posname_list.append(posname)
			last_posname = posname
		PSS[posname] = PSS.get(posname, 0.0) + abs(feature_map[feature])
	with open(str(model_filename).replace("_hypothesis_out_feature_weights", "_PSS"), 'w') as file:
		if len(pos_stats) > 1:
			file.write("{}\t{}\t{}\n".format("Position Name", "PSS", '\t'.join(pos_stats["Position Name"])))
			for posname in posname_list:
				file.write("{}\t{}\t{}\n".format(posname, PSS[posname], '\t'.join(pos_stats[posname])))
		else:
			file.write("{}\t{}\n".format("Position Name", "PSS"))
			for posname in posname_list:
				file.write("{}\t{}\n".format(posname, PSS[posname]))
	# Return sum of all position significance scores
	return sum(PSS.values())


def generate_mapped_weights_file(weights_filename, feature_map_filename, model):
	# Read weights and feature mapping files
	PSS = {}
	posname_list = []
	last_posname = ""
	#model = xml_model_to_dict(weights_filename)
	feature_map = {}
	pos_stats = {}
	# output_filename = str(weights_filename).replace("_hypothesis_out_feature_weights.xml", "_mapped_feature_weights.txt")
	output_filename = str(weights_filename).replace("_out_feature_weights", "_mapped_feature_weights").replace(".xml", ".txt").replace("_hypothesis", "")
	with open(feature_map_filename, 'r') as file:
		for line in file:
			data = line.strip().split("\t")
			if len(data) == 2:
				feature_map[int(data[0])] = data[1]
	if os.path.exists(feature_map_filename.replace("feature_mapping_", "pos_stats_")):
		with open(feature_map_filename.replace("feature_mapping_", "pos_stats_"), 'r') as file:
			for line in file:
				data = line.strip().split("\t")
				if len(data) > 1:
					pos_stats[data[0]] = data[1:]
	with open(output_filename, 'w') as file:
		for i in range(0, len(model["weight_list"])):
			if float(model["weight_list"][i]) != 0.0:
				file.write("{}\t{}\n".format(feature_map[i+1], model["weight_list"][i]))
			posname = feature_map[i + 1][0:-2]
			if posname != last_posname:
				posname_list.append(posname)
				last_posname = posname
			PSS[posname] = PSS.get(posname, 0.0) + abs(model["weight_list"][i])
		file.write("{}\t{}\n".format("Intercept", model["intercept"]))
	with open(str(output_filename).replace("_mapped_feature_weights", "_PSS"), 'w') as file:
		if len(pos_stats) > 1:
			file.write("{}\t{}\t{}\n".format("Position Name", "PSS", '\t'.join(pos_stats["Position Name"])))
			for posname in posname_list:
				file.write("{}\t{}\t{}\n".format(posname, PSS[posname], '\t'.join(pos_stats[posname])))
		else:
			file.write("{}\t{}\n".format("Position Name", "PSS"))
			for posname in posname_list:
				file.write("{}\t{}\n".format(posname, PSS[posname]))
	# Return sum of all position significance scores
	return sum(PSS.values())

