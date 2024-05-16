import os
import sys
import time
import argparse
import shutil
import math
import numpy
import pipeline_funcs as pf
import gene_contribution_visualizer as gcv


# Reads model file into a dictionary
def read_model(filename):
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


# Extracts weights from a list of alignment files based on model dictionary and generate sums per species/gene
def extract_gene_sums(aln_list, model):
	gene_files = {}
	gene_sums = {}
	with open(aln_list, 'r') as file:
		for line in file:
			for aln_filename in line.strip().split(","):
				gene_files[os.path.splitext(os.path.basename(aln_filename.strip()))[0]] = aln_filename.strip()
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
		alignment = read_fasta(gene_files[gene])
		for seq_id in alignment.keys():
			gene_sums[gene][seq_id] = 0
			for pos in model[gene].keys():
				gene_sums[gene][seq_id] += model[gene][pos].get(alignment[seq_id][pos], 0)
	return gene_sums


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


def parse_response_file(response_filename, species_list):
	responses = {seq_id: None for seq_id in species_list}
	with open(response_filename, 'r') as file:
		custom_responses = [tuple(line.strip().split("\t")) for line in file]
		for custom_response in custom_responses:
			if responses[custom_response[0]] is None:
				responses[custom_response[0]] = custom_response[1]
			else:
				raise Exception("Response value of sequence {} specified more than once".format(custom_response[0]))
	# for seq_id in responses.keys():
	# 	if responses[seq_id] is None:
	# 		responses[seq_id] = "0"
	return responses

def main(args):
	model_dir = os.path.dirname(args.model)
	model_basename = os.path.splitext(os.path.basename(args.model))[0].replace("_hypothesis", "").replace("_mapped_feature_weights", "")
	if args.output is None:
		args.output = "{}_applied_gene_prediction.txt".format(model_basename)
	model = read_model(os.path.join(model_dir, model_basename + "_mapped_feature_weights.txt"))
	gene_sums = extract_gene_sums(args.aln_list, model)

	species_list = set()

	for gene in gene_sums.keys():
		species_list.update(list(gene_sums[gene].keys()))

	if args.response is not None:
		response = parse_response_file(args.response, species_list)
	else:
		response = {seq_id: 0 for seq_id in species_list}

	weighted_group_sums = apply_group_weights(args.aln_list, gene_sums, os.path.join(model_dir, "group_indices_" + model_basename + "_hypothesis.txt"), species_list)

	gene_list = list(gene_sums.keys())
	group_list = list(weighted_group_sums.keys())



	with open(args.output, 'w') as file:
		#file.write("SeqID\tResponse\tPrediction\tIntercept\t{}\n".format("\t".join([",".join(x) for x in gene_list])))
		#file.write("SeqID\tResponse\tPrediction\tIntercept\t{}\n".format("\t".join(gene_list)))
		file.write("SeqID\tResponse\tPrediction\tIntercept\t{}\n".format("\t".join(group_list)))


		for seq_id in species_list:
			if response[seq_id] is None:
				continue
			prediction = model["Intercept"]
			for group in group_list:
				prediction += weighted_group_sums[group].get(seq_id, 0)
			file.write("{}\t{}\t{}\t{}\t{}\n".format(seq_id, response[seq_id], prediction, model["Intercept"], "\t".join([str(weighted_group_sums[group].get(seq_id, numpy.nan)) for group in group_list])))

	gcv_files = gcv.main(args.output, gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff, m_grid=args.m_grid)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="ESL Model Application.")
	# should point to a model/hypothesis file inside an ESL results directory, and derive other model file locations from that
	parser.add_argument("model", help="ESL mapped_weights or hypothesis file *inside* an ESL results directory.", type=str)
	# aln_list is part of the novel input set
	parser.add_argument("aln_list", help="List of alignment files to extract features from for labeled positions/alleles.", type=str)
	# response values can be specified, but don't have to be.
	parser.add_argument("--response", help="File containing list of named node/response value pairs, primarily used for sorting visual output.", type=str, default=None)
	parser.add_argument("-o", "--output", help="Output directory.", type=str, default=None)
	parser.add_argument("--gene_display_limit", help="Limits the number of genes displayed in the generated graph images.", type=int, default=100)
	parser.add_argument("--gene_display_cutoff", help="Limits genes displayed in the generated graph images to those with sum-of-squares greater than cutoff value.", type=int, default=0.0)
	parser.add_argument("--m_grid", help="Generate m-grid graphical output.", action='store_true', default=False)
	args = parser.parse_args()
	main(args)










