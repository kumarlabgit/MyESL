import os
import sys
import subprocess
import time
import shutil
import statistics
import random
import psutil
import numpy
import pandas
from datetime import datetime
from datetime import timedelta
from Bio import Phylo
from Bio import AlignIO
import gene_contribution_visualizer as gcv


def generate_hypothesis_set(args):
	if not os.path.exists(args.output):
		os.mkdir(args.output)
	new_files = {'hypothesis_files': [], 'xval_id_files': [], 'slep_opts_files': [], 'sweights_files': []}
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
	responses = {}
	if response_filename is None:
		tree = Phylo.parse(newick_filename, 'newick').__next__()
		tree_tbl = tree.total_branch_length()
		if smart_sampling is not None and tree_tbl == 0:
			smart_sampling = 2
			print("Warning: provided tree does not contain branch lengths, smart sampling will select only/entire sister clade(s) as negative set(s).")
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
				if tree_tbl > 0 and response_sum <= -1:
					negative_set = [key for key in responses[nodename].keys() if responses[nodename][key] == -1]
					for i in range(response_sum, 0):
						# get minimum pair distance from matrix
						# print(negative_set)
						# print(temp_distance.keys())
						# minimum_distance = min([min([temp_distance[t1][t2] for t1 in negative_set if t1 != t2]) for t2 in negative_set])
						min_length = min([nodes[negative_taxa].branch_length for negative_taxa in negative_set])
						min_taxa = set()
						# [[min_taxa.add(t1) for t1 in negative_set if temp_distance[t1][t2] == minimum_distance] for t2 in negative_set]
						[min_taxa.add(t1) for t1 in negative_set if nodes[t1].branch_length == min_length]
						min_taxa = list(min_taxa)
						# Delete the species in min_taxa that occurs first on the tree, assuming tree.get_terminals() is deterministically ordered
						for taxa in taxa_list:
							if taxa in min_taxa:
								responses[nodename][taxa] = 0
								negative_set.remove(taxa)
								break
			elif smart_sampling == 2:
				# print(nodename)
				responses[nodename] = {x: 0 for x in taxa_list}
				for terminal in nodes[nodename].get_terminals():
					responses[nodename][terminal.name] = 1
				target = nodes[nodename]
				response_sum = sum(responses[nodename].values())
				# while response_sum > 0.1 * len(nodes[nodename].get_terminals()):
				while response_sum > 0:
					try:
						parent = tree.get_path(target)[-2]
						for cousin in parent.get_terminals():
							if responses[nodename][cousin.name] == 0:
								responses[nodename][cousin.name] = -1
						response_sum = sum(responses[nodename].values())
						target = parent
						if response_sum < 0.1 * len(nodes[nodename].get_terminals()):
							try:
								parent = tree.get_path(target)[-2]
							except:
								parent = tree.root
							next_sum = 0
							for cousin in parent.get_terminals():
								if responses[nodename][cousin.name] != 1:
									next_sum -= 1
								else:
									next_sum += 1
							if abs(next_sum) > response_sum * 4:
								response_sum = 0
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
				if tree_tbl > 0 and response_sum < 0:
					negative_set = [key for key in responses[nodename].keys() if responses[nodename][key] == -1]
					for i in range(response_sum, 0):
						# get minimum pair distance from matrix
						# print(negative_set)
						# print(temp_distance.keys())
						# minimum_distance = min([min([temp_distance[t1][t2] for t1 in negative_set if t1 != t2]) for t2 in negative_set])
						min_length = min([nodes[negative_taxa].branch_length for negative_taxa in negative_set])
						min_taxa = set()
						# [[min_taxa.add(t1) for t1 in negative_set if temp_distance[t1][t2] == minimum_distance] for t2 in negative_set]
						[min_taxa.add(t1) for t1 in negative_set if nodes[t1].branch_length == min_length]
						min_taxa = list(min_taxa)
						# Delete the species in min_taxa that occurs first on the tree, assuming tree.get_terminals() is deterministically ordered
						for taxa in taxa_list:
							if taxa in min_taxa:
								responses[nodename][taxa] = 0
								negative_set.remove(taxa)
								break
				response_sum = sum(responses[nodename].values())
				#				elif response_sum > 0.1 * len(nodes[nodename].get_terminals()):
				if tree_tbl > 0 and response_sum > 0:
					positive_set = [key for key in responses[nodename].keys() if responses[nodename][key] == 1]
					for i in range(0, response_sum):
						# get minimum pair distance from matrix
						# print(negative_set)
						# print(temp_distance.keys())
						# minimum_distance = min([min([temp_distance[t1][t2] for t1 in positive_set if t1 != t2]) for t2 in positive_set])
						min_length = min([nodes[positive_taxa].branch_length for positive_taxa in positive_set])
						min_taxa = set()
						# [[min_taxa.add(t1) for t1 in positive_set if temp_distance[t1][t2] == minimum_distance] for t2 in positive_set]
						[min_taxa.add(t1) for t1 in positive_set if nodes[t1].branch_length == min_length]
						min_taxa = list(min_taxa)
						# Delete the species in min_taxa that occurs first on the tree, assuming tree.get_terminals() is deterministically ordered
						for taxa in taxa_list:
							if taxa in min_taxa:
								responses[nodename][taxa] = 0
								positive_set.remove(taxa)
								break
			elif smart_sampling == 3:
				# print(nodename)
				responses[nodename] = {x: 0 for x in taxa_list}
				for terminal in nodes[nodename].get_terminals():
					responses[nodename][terminal.name] = 1
				target = nodes[nodename]
				response_sum = sum(responses[nodename].values())
				while response_sum > 0:
					try:
						parent = tree.get_path(target)[-2]
						for cousin in parent.get_terminals():
							if responses[nodename][cousin.name] == 0:
								responses[nodename][cousin.name] = -1
						response_sum = sum(responses[nodename].values())
						target = parent
						if response_sum < 0.1 * len(nodes[nodename].get_terminals()):
							try:
								parent = tree.get_path(target)[-2]
							except:
								parent = tree.root
							next_sum = 0
							for cousin in parent.get_terminals():
								if responses[nodename][cousin.name] != 1:
									next_sum -= 1
								else:
									next_sum += 1
							if abs(next_sum) > response_sum * 4:
								response_sum = 0
					except:
						parent = tree.root
						for cousin in parent.get_terminals():
							if responses[nodename][cousin.name] == 0:
								responses[nodename][cousin.name] = -1
						response_sum = 0
					if len(taxa_list) < 2.0 * len(nodes[nodename].get_terminals()):
						pass
				response_sum = sum(responses[nodename].values())
				if response_sum <= -1:
					negative_set = [key for key in responses[nodename].keys() if responses[nodename][key] == -1]
					for i in range(response_sum, 0):
						# get minimum pair distance from matrix
						# print(negative_set)
						# print(temp_distance.keys())
						# minimum_distance = min([min([temp_distance[t1][t2] for t1 in negative_set if t1 != t2]) for t2 in negative_set])
						min_length = min([nodes[negative_taxa].branch_length for negative_taxa in negative_set])
						min_taxa = set()
						# [[min_taxa.add(t1) for t1 in negative_set if temp_distance[t1][t2] == minimum_distance] for t2 in negative_set]
						[min_taxa.add(t1) for t1 in negative_set if nodes[t1].branch_length == min_length]
						min_taxa = list(min_taxa)
						# Delete the species in min_taxa that occurs first on the tree, assuming tree.get_terminals() is deterministically ordered
						for taxa in taxa_list:
							if taxa in min_taxa:
								responses[nodename][taxa] = 0
								negative_set.remove(taxa)
								break
	else:
		with open(response_filename, 'r') as file:
			basename = os.path.splitext(os.path.basename(response_filename))[0]
			responses[basename] = {}
			taxa_list = []
			custom_responses = [tuple(line.strip().split("\t")) for line in file]
			for response in custom_responses:
				if response[0] in taxa_list:
					raise Exception("Response value of sequence {} specified more than once".format(response[0]))
				taxa_list = taxa_list + [response[0]]
				responses[basename][response[0]] = response[1]
	for nodename in responses.keys():
		pos_idx = 1
		neg_idx = 1
		pos_idxs = []
		neg_idxs = []
		new_fname = os.path.join(args.output, "{}_hypothesis.txt".format(nodename, args.output))
		with open(new_fname, 'w') as file:
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
		new_files['hypothesis_files'] += [new_fname]
		random.shuffle(pos_idxs)
		random.shuffle(neg_idxs)
		new_fname = os.path.join(args.output, "{}_xval_groups.txt".format(nodename, args.output))
		if args.kfold_ids is None:
			with open(new_fname, 'w') as file:
				for taxa in taxa_list:
					if responses[nodename][taxa] not in [0, "0"]:
						if float(responses[nodename][taxa]) > 0:
							# file.write("{}\t{}\n".format(taxa, pos_idxs.pop()))
							file.write("{}\n".format(pos_idxs.pop()))
						elif float(responses[nodename][taxa]) < 0:
							# file.write("{}\t{}\n".format(taxa, neg_idxs.pop()))
							file.write("{}\n".format(neg_idxs.pop()))
		else:
			shutil.copy(args.kfold_ids, new_fname)
		new_files['xval_id_files'] += [new_fname]
		opts_fname = os.path.join(args.output, "{}_slep_opts.txt".format(nodename, args.output))
		if slep_sample_balance:
			sweights_fname = os.path.join(args.output, "{}_sweights.txt".format(nodename, args.output))
			with open(opts_fname, 'w') as opts_file:
				if args.slep_opts is not None:
					with open(args.slep_opts, 'r') as base_opts_file:
						for line in base_opts_file:
							opts_file.write(line)
				# ratio = sum([1.0 for x in responses[nodename].values() if x == 1])/sum([1.0 for x in responses[nodename].values() if x == -1])
				opts_file.write("sWeight\t{}\n".format(os.path.join(args.output, "{}_sweights.txt".format(nodename, args.output))))
				with open(sweights_fname, 'w') as sweights_file:
					sweights_file.write("{}\n".format(1.0))
					sweights_file.write("{}\n".format(
						sum([1.0 for x in responses[nodename].values() if int(x) == 1]) / sum([1.0 for x in responses[nodename].values() if int(x) == -1])))
			new_files['slep_opts_files'] += [opts_fname]
			new_files['sweights_files'] += [sweights_fname]
		elif args.slep_opts is not None:
			shutil.copy(args.slep_opts, opts_fname)
			new_files['slep_opts_files'] += [opts_fname]
	# for file_type in new_files.keys():
	# 	files = []
	# 	for file in new_files[file_type]:
	# 		try:
	# 			shutil.move(file, os.path.join(args.output, file))
	# 			files += [os.path.join(args.output, file)]
	# 		except:
	# 			pass
	# 	new_files[file_type] = files
	if len(new_files['slep_opts_files']) == 0:
		new_files['slep_opts_files'] = [None for i in range(0, len(new_files["hypothesis_files"]))]
	return new_files


def generate_input_matrices(args, file_dict):
	new_files = {'response_files': [], 'group_indices_files': [], 'features_files': [], 'field_files': [], 'feature_mapping_files': [], 'pos_stats_files': []}
	aln_file_list = {}
	gene_list = []
	group_list = []
	partitions_max = 1.0
	output_basename = os.path.split(args.output)[1]
	options = "useCaching threads {}".format(args.threads)
	if not args.include_singletons:
		options = "{} {}".format(options.strip(), "is")  # Note: "is" stands for /ignore/ singletons in the preprocessor.
	if args.upsample_balance:
		options = "{} {}".format(options.strip(), "ub")
	elif args.downsample_balance:
		options = "{} {}".format(options.strip(), "db")
	if args.fuzz_indels:
		options = "{} {}".format(options.strip(), "fuzzIndels")
	if args.bit_ct > 1:
		options = "{} {} {}".format(options.strip(), "ct", args.bit_ct)
	if args.data_type != "universal":
		if args.data_type not in ["nucleotide", "protein", "molecular", "numeric"]:
			raise Exception("Invalid data_type specified ({}), accepted values are: universal, molecular, protein, nucleotide, numeric.".format(data_type))
		options = "{} {} {}".format(options.strip(), "dataType", args.data_type)
	options = options.replace(" ", "*")
	alnlist_dir, alnlist_filename = os.path.split(args.aln_list)
	if alnlist_dir == "":
		alnlist_dir = "."
	preprocess_cwd = os.path.split(args.output)[0]
	if preprocess_cwd == "":
		preprocess_cwd = "."
	if os.name == "posix":
		preprocess_exe = os.path.join(os.getcwd(), "bin", "preprocess")
	else:
		preprocess_exe = os.path.join(os.getcwd(), "bin", "preprocess.exe")
	#for filename in hypothesis_filename_list:
	for filename in file_dict['hypothesis_files']:
		# Construct preprocessing command
		# preprocess_cmd = "{}*{}*{}*{}*{}".format(preprocess_exe, os.path.join(os.getcwd(), filename), alnlist_filename, output_basename, options)
		# preprocess_cmd = "{}*{}*{}*{}*{}".format(preprocess_exe, os.path.join(os.getcwd(), filename), os.path.join(os.getcwd(), args.output, "aln_list.txt"), output_basename, options)
		preprocess_cmd = "{}*{}*{}*{}*{}".format(preprocess_exe, os.path.join(os.getcwd(), filename), os.path.join(os.getcwd(), args.aln_list), output_basename, options)
		print(preprocess_cmd.replace("*"," "))
		hypothesis_basename = os.path.splitext(os.path.basename(filename))[0]
		if args.skip_preprocessing:
			if os.path.exists(os.path.join(args.output, "feature_" + hypothesis_basename + ".txt")):
				print("Features file detected, skipping preprocessing step...")
			else:
				raise Exception("Preprocessing skipped, but no features file detected at {}.".format(os.path.join(args.output, "feature_" + hypothesis_basename + ".txt")))
		else:
			subprocess.call(preprocess_cmd.split("*"), stderr=subprocess.STDOUT, cwd=preprocess_cwd)
			if not args.disable_mc:
				partitions_max = max(partitions_max, check_memory("feature_{}.txt".format(hypothesis_basename), args.output))
			if args.data_type != "numeric" and partitions_max <= 1.0:
				print("Calculating position statistics...")
				position_stats = {}
				#stat_keys = ["mic", "entropy"]
				stat_keys = ["mic"]
				pos_count = 0
				for aln_basename in aln_file_list.keys():
					position_stats[aln_basename] = calculate_position_stats(aln_file_list[aln_basename], filename)
					pos_count = pos_count + 1
					if pos_count % 1000 == 0:
						print("Calculated position stats for {}/{} input files...".format(pos_count, len(aln_file_list)))
				with open(os.path.join(args.output, "pos_stats_" + hypothesis_basename + ".txt"), 'w') as file:
					file.write("{}\t{}\n".format("Position Name", '\t'.join(stat_keys)))
					for aln_basename in position_stats.keys():
						for i in range(0, len(position_stats[aln_basename]["mic"])):
							file.write("{}_{}\t{}\n".format(aln_basename, i, '\t'.join([str(position_stats[aln_basename][stat_key][i]) for stat_key in stat_keys])))
				new_files['pos_stats_files'] += [os.path.join(args.output, "pos_stats_" + hypothesis_basename + ".txt")]
			shutil.move(os.path.join(args.output, "feature_" + output_basename + ".txt"), os.path.join(args.output, "feature_" + hypothesis_basename + ".txt"))
			shutil.move(os.path.join(args.output, "group_indices_" + output_basename + ".txt"), os.path.join(args.output, "group_indices_" + hypothesis_basename + ".txt"))
			shutil.move(os.path.join(args.output, "response_" + output_basename + ".txt"), os.path.join(args.output, "response_" + hypothesis_basename + ".txt"))
#			shutil.move(hypothesis_basename.replace("hypothesis", "xval_groups.txt"), os.path.join(oargs.output, hypothesis_basename.replace("hypothesis", "xval_groups.txt")))
			shutil.move(os.path.join(args.output, "field_" + output_basename + ".txt"), os.path.join(args.output, "field_" + hypothesis_basename + ".txt"))
			shutil.move(os.path.join(args.output, "feature_mapping_" + output_basename + ".txt"), os.path.join(args.output, "feature_mapping_" + hypothesis_basename + ".txt"))
#			try:
#				shutil.move(os.path.join(args.output, "resampled_" + output_basename + ".txt"), os.path.join(args.output, "resampled_" + hypothesis_basename + ".txt"))
#			except:
#				pass
		#new_files = {'response_files': [], 'group_indices_files': [], 'features_files': [], 'field_files': []}
		#response_file_list.append(os.path.join(args.output, "response_" + hypothesis_basename + ".txt"))
		new_files['response_files'] += [os.path.join(args.output, "response_" + hypothesis_basename + ".txt")]
		#group_indices_file_list.append(os.path.join(args.output, "group_indices_" + hypothesis_basename + ".txt"))
		new_files['group_indices_files'] += [os.path.join(args.output, "group_indices_" + hypothesis_basename + ".txt")]
		#field_file_list.append(os.path.join(args.output, "field_" + hypothesis_basename + ".txt"))
		new_files['field_files'] += [os.path.join(args.output, "field_" + hypothesis_basename + ".txt")]
		#features_file_list.append(os.path.join(args.output, "feature_" + hypothesis_basename + ".txt"))
		new_files['features_files'] += [os.path.join(args.output, "feature_" + hypothesis_basename + ".txt")]
		new_files['feature_mapping_files'] += [os.path.join(args.output, "feature_mapping_" + hypothesis_basename + ".txt")]
	if partitions_max > 2.0:
		raise Exception("Minimum required subsets:{}".format(partitions_max))
	with open(os.path.join(args.output, "missing_seqs_" + output_basename + ".txt"), "r") as file:
		for ms_line in file:
			ms_data = ms_line.strip().split("\t")
			args.missing_seqs.add((ms_data[1], os.path.splitext(os.path.basename(ms_data[0]))[0]))
	return new_files


def generate_lambda_list(args, file_dict):
	new_files = {}
	lambda_list = []
	[z_min, z_max, z_interval] = [float(x) for x in args.grid_z.strip().split(',')]
	z_steps = (z_max - z_min) / z_interval
	z_list = [z_min + (x * z_interval) for x in range(0, int(z_steps) + 1) if z_min + (x * z_interval) < 1.0]
	[y_min, y_max, y_interval] = [float(x) for x in args.grid_y.split(',')]
	y_steps = (y_max - y_min) / y_interval
	y_list = [y_min + (x * y_interval) for x in range(0, int(y_steps) + 1) if y_min + (x * y_interval) < 1.0]
	with open(os.path.join(args.output, "lambda_list.txt".format(args.output)), 'w') as file:
		for lambda1 in z_list:
			for lambda2 in y_list:
				lambda_list.append((lambda1, lambda2))
				file.write("{:g}\t{:g}\n".format(lambda1, lambda2))
	new_files['lambda_list_file'] = [os.path.join(args.output, "lambda_list.txt".format(args.output))]
	return new_files


def run_sg_lasso(args, file_dict):
	new_files = {"weights_files": []}
	if args.method == "leastr":
		method = "sg_lasso_leastr"
	elif args.method == "logistic":
		method = "sg_lasso"
	elif args.method == "ol_leastr":
		method = "overlapping_sg_lasso_leastr"
	elif args.method == "ol_logistic":
		method = "overlapping_sg_lasso_logisticr"
	else:
		raise Exception("Provided method name not recognized, please provide a valid method name.")
	if os.name == "posix":
		esl_exe = os.path.join(os.getcwd(), "bin", "{}".format(method))
	else:
		esl_exe = os.path.join(os.getcwd(), "bin", "{}.exe".format(method))
	lambda_list = []
	with open(file_dict["lambda_list_file"][0], 'r') as file:
		for line in file:
			data = line.strip().split("\t")
			lambda_list.append((data[0], data[1]))
	for response_filename, features_filename, groups_filename, field_filename, slep_opts_filename in \
			zip(file_dict["response_files"], file_dict["features_files"], file_dict["group_indices_files"], file_dict["field_files"], file_dict["slep_opts_files"]):
		basename = str(os.path.splitext(os.path.basename(response_filename))[0]).replace("response_", "")
		if slep_opts_filename is None:
			esl_cmd = "{}*-f*{}*-z*{}*-y*{}*-n*{}*-r*{}*-w*{}".format(esl_exe, features_filename, args.lambda1, args.lambda2, groups_filename, response_filename, os.path.join(args.output, basename + "_out_feature_weights"))
		else:
			esl_cmd = "{}*-f*{}*-z*{}*-y*{}*-n*{}*-r*{}*-s*{}*-w*{}".format(esl_exe, features_filename, args.lambda1, args.lambda2, groups_filename, response_filename, slep_opts_filename, os.path.join(args.output, basename + "_out_feature_weights"))
		if method in ["overlapping_sg_lasso_leastr", "overlapping_sg_lasso_logisticr"]:
			esl_cmd = esl_cmd + "*-g*{}".format(field_filename)
		esl_cmd = esl_cmd + "*-l*{}".format(file_dict["lambda_list_file"][0])
		if args.grid_gene_threshold > 0:
			esl_cmd = esl_cmd + "*-c*{}".format(args.grid_gene_threshold)
		esl_cmd = esl_cmd + "*--model_format*flat"
		print(esl_cmd.replace("*", " "))
		subprocess.call(esl_cmd.split("*"), stderr=subprocess.STDOUT)
		new_files["weights_files"].append([os.path.join(args.output, "{}_out_feature_weights_{}_{}.txt".format(basename, lambda_val2label(val[0]), lambda_val2label(val[1]))) for val in lambda_list])
	return new_files


def process_weights(args, file_dict):
	new_files = {"gene_prediction_files": [], "HSS_file": [], "GSS_files": [], "PSS_files": []}
	aln_lib = {}
	for (weights_filename_list, hypothesis_filename, groups_filename) in zip(file_dict["weights_files"], file_dict["hypothesis_files"], file_dict["group_indices_files"]):
		missing_results = []
		gene_prediction_files = []
		for weights_filename in weights_filename_list:
			outname = weights_filename.replace("_hypothesis", "").replace("_out_feature_weights", "_gene_predictions").replace(".xml", ".txt")
			if os.path.exists(weights_filename):
				args.HSS[os.path.basename(hypothesis_filename)] = args.HSS.get(os.path.basename(hypothesis_filename), 0) + process_single_grid_weight(weights_filename, hypothesis_filename, groups_filename, outname, aln_lib, args)
				gene_prediction_files += [outname]
			elif args.grid_gene_threshold > 0:
				print("No results file detected, most likely due to grid_gene_threshold: {}".format(weights_filename))
				missing_results.append(weights_filename)
			else:
				raise Exception("Missing results file: {}".format(weights_filename))
		for missing_result in missing_results:
			weights_filename_list.remove(missing_result)
		new_files["gene_prediction_files"].append(gene_prediction_files)
		new_files["GSS_files"].append([fname.replace("gene_predictions", "GSS") for fname in gene_prediction_files])
		new_files["PSS_files"].append([fname.replace("gene_predictions", "PSS") for fname in gene_prediction_files])
	with open(os.path.join(args.output, "HSS.txt"), 'w') as file:
		file.write("{}\t{}\n".format("Hypothesis", "HSS"))
		for hypothesis_filename in [os.path.basename(filename) for filename in file_dict["hypothesis_files"]]:
			file.write("{}\t{}\n".format(hypothesis_filename.replace("_hypothesis.txt", ""), args.HSS[os.path.basename(hypothesis_filename)]))
	new_files["HSS_file"] += [os.path.join(args.output, "HSS.txt")]
	return new_files


def process_single_grid_weight(weights_filename, hypothesis_filename, groups_filename, outname, aln_lib, args):
	numeric = False
	if args.data_type == "numeric":
		numeric = True
	apply_ESL_model(args.aln_list, aln_lib, weights_filename, hypothesis_filename, groups_filename, outname, args.missing_seqs, numeric=numeric)
	total_significance = generate_significance_scores(weights_filename, groups_filename, numeric=numeric)
	return total_significance


def apply_ESL_model(aln_list, aln_lib, model_file, hypothesis_file, groups_filename, output_filename, missing_seqs, numeric=False):
	model = read_ESL_model(model_file, numeric=numeric)
	gene_sums, gene_significance_scores = extract_gene_sums(aln_list, aln_lib, model, numeric=numeric)
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
		file.write("SeqID\tResponse\tPrediction\tIntercept\t{}\n".format("\t".join(group_list)))
		for seq_id in list(species_list) + missing_seqids:
			if response[seq_id] is None:
				continue
			prediction = model["Intercept"]
			for group in group_list:
				prediction += weighted_group_sums[group].get(seq_id, 0)
			for group in group_list:
				if sum([1 for x in group.split(",") if (seq_id, x) in missing_seqs]) == len(group.split(",")):
					weighted_group_sums[group][seq_id] = numpy.nan
			file.write("{}\t{}\t{}\t{}\t{}\n".format(seq_id, response[seq_id], prediction, model["Intercept"],
													 "\t".join([str(weighted_group_sums[group].get(seq_id, numpy.nan)) for group in group_list])))
	if "gene_predictions_xval" not in output_filename:
		with open(str(output_filename).replace("_gene_predictions", "_GSS"), 'w') as file:
			file.write("{}\t{}\n".format("Gene", "GSS"))
			for gene in gene_significance_scores.keys():
				file.write("{}\t{}\n".format(gene, str(gene_significance_scores[gene])))


def read_ESL_model(filename, numeric=False):
	model = {}
	last_gene = ""
	last_pos = -1
	if numeric:
		with open(filename, 'r') as file:
			for line in file:
				data = line.strip().split("\t")
				if data[0] == "Intercept":
					model["Intercept"] = float(data[1])
					continue
				feature = data[0].split("_")
				gene = "_".join(feature[0:-1])
				pos = int(feature[-1])
				weight = float(data[1])
				if gene != last_gene:
					model[gene] = {pos: weight}
				else:
					model[gene].update({pos: weight})
				last_gene = gene
		return model
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


def extract_gene_sums(aln_list, aln_lib, model, numeric=False):
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
			if numeric:
				aln_lib[gene] = read_numeric(gene_files[gene])
			else:
				aln_lib[gene] = read_fasta(gene_files[gene])
		for seq_id in aln_lib[gene].keys():
			if numeric:
				gene_sums[gene][seq_id] = sum([model[gene][pos] * aln_lib[gene][seq_id][pos] for pos in model[gene].keys()])
			else:
				gene_sums[gene][seq_id] = sum([model[gene][pos].get(aln_lib[gene][seq_id][pos], 0) for pos in model[gene].keys()])
		# gene_signifcance_scores[gene] = sum([sum(model[gene][pos].values()) for pos in model[gene].keys()])
		if numeric:
			gene_signifcance_scores[gene] = sum([abs(val) for val in model[gene].values()])
		else:
			gene_signifcance_scores[gene] = sum([sum([abs(val) for val in pos.values()]) for pos in model[gene].values()])
	return gene_sums, gene_signifcance_scores


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


def read_numeric(filename, delimiter='\t', header=False):
	with open(filename, 'r') as file:
		values = {}
		for line in file:
			if header:
				header = False
				continue
			data = line.strip().split(delimiter)
			values[data[0]] = [float(x) for x in data[1:]]
	return values


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


def generate_significance_scores(model_filename, groups_filename, numeric=False):
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
	if numeric:
		with open(str(model_filename).replace("_hypothesis_out_feature_weights", "_PSS"), 'w') as file:
			file.write("{}\t{}\n".format("Position Name", "PSS"))
			for posname in feature_map.keys():
				file.write("{}\t{}\n".format(posname, abs(feature_map[posname])))
		return sum([abs(x) for x in feature_map.values()])
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


def generate_model_graphics(args, file_dict):
	new_files = {"GCV_files": [], "SPS_SPP_files": []}
	for gene_prediction_file_list in file_dict["gene_prediction_files"]:
		GCV_files = []
		for gene_prediction_file in gene_prediction_file_list:
			try:
				GCV_files += [gcv.main(gene_prediction_file, species_limit=args.species_display_limit, gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff, m_grid=args.m_grid)]
			except:
				print("Could not produce graphic for {}.".format(gene_prediction_file))
		SPS_SPP_files = [fname.replace("gene_predictions", "SPS_SPP").replace(".png", ".txt") for fname in GCV_files]
		new_files["GCV_files"].append(GCV_files)
		new_files["SPS_SPP_files"].append(SPS_SPP_files)
	return new_files


def summarize_models(args, file_dict):
	new_files = {"GCS_median_files": [], "GSS_median_files": [], "PSS_median_files": [], "GCV_median_files": [], "SPS_SPP_median_files": []}
	pss_vals = {}
	gss_vals = {}
	for hypothesis_file, gene_prediction_file_list in zip(file_dict["hypothesis_files"], file_dict["gene_prediction_files"]):
		gene_predictions = []
		for gene_prediction_file in gene_prediction_file_list:
			rmse = 0
			acc = 1
			model = pandas.read_csv(gene_prediction_file, delim_whitespace=True, index_col=0)
			# Calculate acc and RMSE of model
			rmse = ((model.Response - model.Prediction) ** 2).mean() ** 0.5
			acc = sum([1 for x in zip(model.Response, model.Prediction) if (x[0] == -1 and x[1] < 0) or (x[0] == 1 and x[1] > 0)]) / float(model.shape[0])
			if rmse <= args.grid_rmse_cutoff and acc >= args.grid_acc_cutoff:
				gene_predictions.append(model)
			# Parse GSS file into GSS vals
			with open(gene_prediction_file.replace("gene_predictions", "GSS"), 'r') as gss_file:
				for line in gss_file:
					data = line.strip().split('\t')
					if data[0] == "Gene":
						continue
					if data[0] in gss_vals:
						gss_vals[data[0]].update({gss_file: data[1]})
					else:
						gss_vals[data[0]] = {gss_file: data[1]}
			with open(gene_prediction_file.replace("gene_predictions", "PSS"), 'r') as pss_file:
				for line in pss_file:
					data = line.strip().split('\t')
					if data[0] == "Position Name" or float(data[1]) == 0.0:
						continue
					if data[0] in pss_vals:
						pss_vals[data[0]].update({pss_file: data[1]})
					else:
						pss_vals[data[0]] = {pss_file: data[1]}
		# Write aggregated GSS file
		new_files["GSS_median_files"] += [hypothesis_file.replace("_hypothesis.txt", "_GSS_median.txt")]
		with open(hypothesis_file.replace("_hypothesis.txt", "_GSS_median.txt"), 'w') as file:
			for gene in gss_vals.keys():
				gss_temp = [float(x) for x in gss_vals[gene].values() if abs(float(x)) > 0]
				if len(gss_temp) == 0:
					gss_temp = [0]
				file.write("{}\t{}\n".format(gene, statistics.median(gss_temp)))
		# Write aggregated PSS file
		new_files["PSS_median_files"] += [hypothesis_file.replace("_hypothesis.txt", "_PSS_median.txt")]
		with open(hypothesis_file.replace("_hypothesis.txt", "_PSS_median.txt"), 'w') as file:
			for pos in pss_vals.keys():
				pss_temp = [float(x) for x in pss_vals[pos].values() if abs(float(x)) > 0]
				if len(pss_temp) == 0:
					pss_temp = [0]
				file.write("{}\t{}\n".format(pos, statistics.median(pss_temp)))
		# Write aggregated gene_predictions file
		new_files["GCS_median_files"] += [hypothesis_file.replace("_hypothesis.txt", "_GCS_median.txt")]
		with open(hypothesis_file.replace("_hypothesis.txt", "_GCS_median.txt"), 'w') as file:
			gene_list = []
			species_list = set()
			for gene_prediction in gene_predictions:
				species_list.update(gene_prediction.index)
				for column in gene_prediction.columns:
					if column not in gene_list:
						gene_list.append(column)
			file.write("SeqID\tResponse\tPrediction_mean\t{}\n".format("\t".join([gene for gene in gene_list if gene not in ["SeqID", "Prediction", "Intercept", "Response"]])))
			species_list = list(species_list)
			for species in species_list:
				file.write("{}\t".format(species))
				for gene in gene_list:
					if gene in ["SeqID", "Intercept"]:
						continue
					GCS_list = [x[gene][species] for x in gene_predictions if gene in x.columns.values.tolist() and species in x[gene].keys() and x[gene][species] != 0]
					if len(GCS_list) == 0:
						GCS_list = [0]
					if gene == "Prediction":
						#file.write("{}\t".format(statistics.mode(GCS_list)))
						file.write("{}\t".format(statistics.mean(GCS_list)))
					else:
						file.write("{}\t".format(statistics.median(GCS_list)))
				file.write("\n")
		new_files["GCV_median_files"] += [gcv.main(hypothesis_file.replace("_hypothesis.txt", "_GCS_median.txt"), lead_cols=3, species_limit=args.species_display_limit, gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff, m_grid=args.m_grid)]
	new_files["SPS_SPP_median_files"] = [fname.replace("GCS_median", "SPS_SPP_median").replace(".png", ".txt") for fname in new_files["GCV_median_files"]]
	return new_files


def grid_search(args):
	file_dict = {}
	args.HSS = {}
	args.missing_seqs = set()
	try:
		if not os.path.exists(args.output):
			os.mkdir(args.output)
		args.timers["preprocessing"]["start"] = datetime.now()
		new_files = generate_hypothesis_set(args)
		file_dict.update(new_files)
		new_files = generate_lambda_list(args, file_dict)
		file_dict.update(new_files)
		new_files = generate_input_matrices(args, file_dict)
		file_dict.update(new_files)
		args.timers["preprocessing"]["total"] = args.timers["preprocessing"].get("total", timedelta(seconds=0)) + (datetime.now() - args.timers["preprocessing"]["start"])
		args.timers["sglasso"]["start"] = datetime.now()
		new_files = run_sg_lasso(args, file_dict)
		file_dict.update(new_files)
		args.timers["sglasso"]["total"] = args.timers["sglasso"].get("total", timedelta(seconds=0)) + (datetime.now() - args.timers["sglasso"]["start"])
		args.timers["analysis"]["start"] = datetime.now()
		new_files = process_weights(args, file_dict)
		file_dict.update(new_files)
		if not args.grid_summary_only:
			new_files = generate_model_graphics(args, file_dict)
			file_dict.update(new_files)
		new_files = summarize_models(args, file_dict)
		file_dict.update(new_files)
		args.timers["analysis"]["total"] = args.timers["analysis"].get("total", timedelta(seconds=0)) + (datetime.now() - args.timers["analysis"]["start"])
		# for file_type in file_dict.keys():
		# 	print(file_type)
		# 	for file in file_dict[file_type]:
		# 		if file is not None:
		# 			if len(file[0]) == 1:
		# 				print("\t{}".format(file))
		# 			else:
		# 				for nested_file in file:
		# 					print("\t\t{}".format(nested_file))
		check_file_dict(file_dict)
		cleanup_directory(args, file_dict)
	except Exception as e:
		raise Exception(e)
	finally:
		message_triggered = False
		for file_type in file_dict.keys():
			for file in file_dict[file_type]:
				#print(file)
				try:
					pass
					#shutil.move(file, os.path.join(args.output, file))
				except:
					pass
	return file_dict


def lambda_val2label(lambda_val):
	lambda_val = float(lambda_val)
	if "{:g}".format(lambda_val)[0:2] == "0.":
		return "{:g}".format(lambda_val)[2:]
	else:
		return "{:g}".format(lambda_val)


def lookup_by_names(tree):
	names = {}
	for clade in tree.find_clades():
		if clade.name:
			if clade.name in names:
				raise ValueError("Duplicate key: %s" % clade.name)
			names[clade.name] = clade
	return names


def check_file_dict(file_dict):
	for file_type in file_dict.keys():
		for file in file_dict[file_type]:
			if file is not None:
				if len(file[0]) > 1:
					for nested_file in file:
						if not os.path.exists(nested_file):
							print("{} not found.".format(nested_file))
				elif not os.path.exists(file):
					print("{} not found.".format(file))


def check_memory(features_filename, output_name):
	#stats_filename = features_filename.replace("feature_", "feature_stats_")
	stats_filename = os.path.join(output_name, "feature_stats_{}.txt".format(os.path.split(output_name)[1]))
	try:
		stats = {}
		with open(stats_filename, 'r') as stats_file:
			for line in stats_file:
				data = line.strip().split("\t")
				stats[data[0]] = data[1]
	except:
		raise Exception("Problem opening feature stats file {}.".format(stats_filename))
	available_mem = psutil.virtual_memory().available
	feature_mem = round(7 * int(stats["Samples"]) * int(stats["Features"]))
	# feature_mem = round((12*10*5) * int(stats["Samples"]) * int(stats["Features"]))
	if available_mem < feature_mem:
		msg = "Exceeding available memory will severely degrade performance, if you're sure you want to try anyways, rerun MyESL with the --disable_mc flag."
		print("Total size of {} in memory ({} bytes) would exceed available memory of {} bytes.\n{}".format(features_filename, feature_mem, available_mem, msg))
		#raise Exception("Total size of {} in memory ({} bytes) would exceed available memory of {} bytes.\n{}".format(features_filename, feature_mem, available_mem, msg))
	return max(2.0, float(feature_mem * 2)/float(available_mem))


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
		#print(aln_filename)
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


def stat_aln_files(aln_list):
	file_stats = {}
	aln_list_dir = os.path.split(aln_list)[0]
	with open(aln_list) as file:
		for line in file:
			group_size = 0
			for aln_filename in line.strip().split(","):
				group_size += os.path.getsize(os.path.join(aln_list_dir, aln_filename))
			file_stats[line.strip()] = group_size
	return file_stats


def generate_aln_partitions(aln_stats, partitions, aln_list, combo_size=2):
	if partitions == combo_size:
		raise Exception("Partition count({}) must be larger than combo_size({}).".format(partitions, combo_size))
	dir, basename = os.path.split(aln_list)
	group_lists = [[] for i in range(0, partitions)]
	group_sizes = [0 for i in range(0, partitions)]
	aln_list_partition_files = []
	shuffled_groups = list(aln_stats.keys())
	# For the sake of reproducibility, don't actually shuffle
	# random.shuffle(shuffled_groups)
	group_idx = 0
	if len(shuffled_groups) > 20 * partitions:
		while group_idx < 0.7 * len(shuffled_groups):
			for i in range(0, partitions):
				group_lists[i] += [shuffled_groups[group_idx]]
				group_sizes[i] += aln_stats[shuffled_groups[group_idx]]
				group_idx += 1
	while group_idx < len(shuffled_groups):
		min_group_idx = group_sizes.index(min(group_sizes))
		group_lists[min_group_idx] += [shuffled_groups[group_idx]]
		group_sizes[min_group_idx] += aln_stats[shuffled_groups[group_idx]]
		group_idx += 1
	# for i in range(0, partitions):
	combo_count = (partitions / 2) * (partitions - 1)
	combo_id = 0
	for i in range(partitions-1, 0, -1):
		for j in range(i-1, -1, -1):
			aln_list_partition_files += [os.path.join(dir, "{}_part{}{}".format(os.path.splitext(basename)[0], combo_id, os.path.splitext(basename)[1]))]
			with open(aln_list_partition_files[combo_id], 'w') as file:
				for group in group_lists[i]:
					file.write("{}\n".format(group))
				for group in group_lists[j]:
					file.write("{}\n".format(group))
			combo_id += 1
	# for i in range(0, combo_count):
	# 	aln_list_partition_files += [os.path.join(dir, "{}_part{}{}".format(os.path.splitext(basename)[0], i, os.path.splitext(basename)[1]))]
	# 	with open(aln_list_partition_files[i], 'w') as file:
	# 		for group in group_lists[i]:
	# 			file.write("{}\n".format(group))
	return aln_list_partition_files


def select_top_gss(GSS_file, pct_select=80.0):
	GSS_list = []
	with open(GSS_file, 'r') as file:
		for line in file:
			GSS_list.append(line.strip().split('\t'))
	GSS_list = [[x[0], float(x[1])] for x in GSS_list]
	GSS_list.sort(key=lambda tup: tup[1], reverse=True)
	select_cnt = 0
	select_sum = 0
	select_list = []
	total_GSS = sum([x[1] for x in GSS_list])
	while select_sum < total_GSS * (pct_select / 100.0):
		select_list += [GSS_list[select_cnt][0]]
		select_sum += GSS_list[select_cnt][1]
		select_cnt += 1
	return select_list


def cleanup_directory(args, file_dict):
	args.verbose = False
	if args.single_lambda_pair:
		for summary_type in [output_type for output_type in file_dict.keys() if "median" in output_type]:
			for fname in file_dict[summary_type]:
				try:
					os.remove(fname)
				except:
					pass
	for output_type in [mode for mode in ["B", "P", "G", "H"] if mode not in args.stats_out]:
		for fname_list in file_dict.get("{}SS_files".format(output_type), []):
			for fname in fname_list:
				try:
					os.remove(fname)
				except:
					pass
		for fname in file_dict.get("{}SS_median_files".format(output_type), []):
			try:
				os.remove(fname)
			except:
				pass
	if "S" not in args.stats_out:
		for fname in file_dict.get("SPS_SPP_median_files", []):
			try:
				os.remove(fname)
			except:
				pass
		for fname_list in file_dict.get("SPS_SPP_files", []):
			for fname in fname_list:
				try:
					os.remove(fname)
				except:
					pass
	for file_type in ["SPS_SPP", "PSS", "GSS", "hypothesis"]:
		for fname_list in file_dict.get("{}_files".format(file_type), []):
			for fname in fname_list:
				if len(fname_list[0]) == 1:
					fname = fname_list
				dir = os.path.split(fname)[0]
				basename = os.path.split(fname)[1]
				new_basename = "{}_{}".format(file_type.replace("_files", ""), basename.replace("_{}".format(file_type).replace("_files", ""), "").replace("_{}".format(os.path.split(fname)[0]), ""))
				# new_fname = "{}_{}".format(file_type, fname.replace("_{}".format(file_type), "").replace("_{}".format(os.path.split(fname)[0]), ""))
				new_fname = os.path.join(dir, new_basename)
				if args.verbose:
					print("Attempting to move {} to {} (or delete it)".format(fname, new_fname))
				try:
					if args.grid_summary_only and file_type != "hypothesis":
						os.remove(fname)
					else:
						shutil.move(fname, new_fname)
				except:
					pass
				if len(fname_list[0]) == 1:
					break
	for fname in file_dict.get("GCS_median_files", []):
		new_fname = os.path.join(os.path.split(fname)[0], "M-Grid_{}".format(os.path.split(fname)[1].replace("_GCS_median".format(args.output), "")))
		if args.verbose:
			print("Attempting to move {} to {}".format(fname, new_fname))
		try:
			shutil.move(fname, new_fname)
			shutil.move(fname.replace(".txt", ".png"), new_fname.replace(".txt", ".png"))
		except:
			pass
	for file_type in [key for key in file_dict.keys() if "median" in key and "GCS" not in key and "GCV" not in key]:
		new_fnames = []
		for fname in file_dict.get(file_type, []):
			dir = os.path.split(fname)[0]
			basename = os.path.split(fname)[1]
			if basename[0:3] == file_type[0:3]:
				continue
			new_basename = "{}_{}".format(file_type.replace("_median_files", ""), basename.replace("_{}".format(file_type).replace("_files", ""), "_summary").replace("_{}".format(os.path.split(fname)[0]), ""))
			# new_fname = "{}_{}".format(file_type.replace("_median_files", ""), fname.replace(file_type.replace("_files", ""), "summary"))
			new_fname = os.path.join(dir, new_basename)
			new_fnames += [new_fname]
			if args.verbose:
				print("Attempting to move {} to {}".format(fname, new_fname))
			try:
				shutil.move(fname, new_fname)
			except:
				pass
		file_dict[file_type] = new_fnames
	for fname_list in file_dict.get("weights_files", []):
		for fname in fname_list:
			new_fname = os.path.join(os.path.split(fname)[0], "MyESL_model_{}".format(os.path.split(fname)[1].replace("_hypothesis_out_feature_weights".format(os.path.split(fname)[0]), "")))
			if args.verbose:
				print("Attempting to move {} to {} (or delete it).".format(fname, new_fname))
			try:
				if args.grid_summary_only:
					os.remove(fname)
				else:
					shutil.move(fname, new_fname)
			except:
				pass
	for fname_list in file_dict.get("gene_prediction_files", []):
		for fname in fname_list:
			for temp_fname in [fname, "{}.png".format(os.path.splitext(fname)[0])]:
				new_fname = os.path.join(os.path.split(temp_fname)[0], "GSC_{}".format(os.path.split(temp_fname)[1].replace("_gene_predictions".format(os.path.split(temp_fname)[0]), "")))
				if args.verbose:
					print("Attempting to move {} to {} (or delete it)".format(fname, new_fname))
				try:
					if args.grid_summary_only:
						os.remove(fname)
					else:
						shutil.move(temp_fname, new_fname)
				except:
					pass
	if not args.preserve_inputs:
		for output_type in ["features", "feature_mapping", "field", "pos_stats", "response", "xval_id", "group_indices", "slep_opts", "sweights"]:
			for fname in file_dict.get("{}_files".format(output_type), []):
				try:
					os.remove(fname)
				except:
					pass


def aln_list_to_absolute(aln_list, output):
	# Convert all alignment paths in aln_list to absolute paths and write new alignment list to output/aln_list.txt
	new_aln_list_filename = os.path.join(output, "aln_list.txt")
	aln_file_list = {}
	alnlist_dir, alnlist_filename = os.path.split(aln_list)
	if alnlist_dir == "":
		alnlist_dir = "."
	alnlist_dir_abspath = os.path.abspath(alnlist_dir)
	with open(new_aln_list_filename, 'w') as abspath_aln_list:
		with open(aln_list) as file:
			for line in file:
				abspath_group = []
				for aln_filename in line.strip().split(","):
					basename = os.path.splitext(os.path.basename(aln_filename.strip()))[0]
					if os.name == "posix":
						aln_filename = aln_filename.replace("\\", "/")
					else:
						aln_filename = aln_filename.replace("/", "\\")
					aln_abspath = os.path.join(alnlist_dir_abspath, os.path.split(aln_filename.strip())[0], os.path.split(aln_filename.strip())[1])
					if basename in aln_file_list.keys():
						if aln_file_list[basename] != aln_abspath:
							raise Exception("Found multiple alignment files with identical basename {}.".format(basename))
					else:
						# aln_file_list[basename] = aln_filename.strip()
						aln_file_list[basename] = aln_abspath
						abspath_group.append(aln_abspath)
				abspath_aln_list.write("{}\n".format(",".join(abspath_group)))
	return new_aln_list_filename


