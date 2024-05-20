import os
import sys
import time
import argparse
import shutil
import math
import statistics
import copy
import pipeline_funcs as pf
import pandas
import gene_contribution_visualizer as gcv
from multiprocessing import Process
from multiprocessing import Pool
from datetime import datetime

timers = {}

def grid_search(args, original_output, input_files):
	hypothesis_file_list = input_files["hypothesis_file_list"]
	slep_opts_file_list = input_files["slep_opts_file_list"]
	features_filename_list = input_files["features_filename_list"]
	groups_filename_list = input_files["groups_filename_list"]
	response_filename_list = input_files["response_filename_list"]
	gene_list = input_files["gene_list"]
	field_filename_list = input_files["field_filename_list"]
	group_list = input_files["group_list"]
	lambda_list_filename = input_files["lambda_list"]
	HSS = {}
	missing_seqs = set()
	gcv_files = []
	# features_filename_list, groups_filename_list, response_filename_list, gene_list, field_filename_list, group_list = pf.generate_input_matrices(args.aln_list, hypothesis_file_list, args)
	with open(os.path.join(original_output, "missing_seqs_" + original_output + ".txt"), "r") as file:
		for ms_line in file:
			ms_data = ms_line.strip().split("\t")
			missing_seqs.add((ms_data[1], os.path.splitext(os.path.basename(ms_data[0]))[0]))
	start = datetime.now()
	weights_file_list = pf.run_esl_grid(features_filename_list, groups_filename_list, response_filename_list, field_filename_list, args.lambda1, args.lambda2, args.method, slep_opts_file_list, lambda_list_filename, args)
	timers["LASSO"] = datetime.now() - start
	print("Time elapsed while performing LASSO operations: {}".format(timers["LASSO"]))
	start = datetime.now()
	missing_results = pf.process_sparse_grid_weights(weights_file_list, hypothesis_file_list, groups_filename_list, HSS, missing_seqs, args)
	timers["post_processing"] = datetime.now() - start
	print("Total time elapsed while processing LASSO results: {}".format(timers["post_processing"]))
	#os.mkdir(args.output)
	lambda_list = []
	missing_predictions = []
	with open(lambda_list_filename, 'r') as file:
		for lambda_line in file:
			lambda_pair1 = lambda_line.strip().split("\t")
			lambda_list.append((str(lambda_pair1[0]).replace("0.", ""), str(lambda_pair1[1]).replace("0.", "")))
	for hypothesis_filename in hypothesis_file_list:
		for lambda_pair1 in lambda_list:
			if hypothesis_filename.replace(".txt","_out_feature_weights_{}_{}".format(lambda_pair1[0], lambda_pair1[1])) in [os.path.splitext(missing_results_fname)[0] for missing_results_fname in missing_results]:
				missing_predictions.append(hypothesis_filename.replace("hypothesis.txt", "gene_predictions_{}_{}.txt".format(lambda_pair1[0], lambda_pair1[1])))
				continue
			# shutil.move(hypothesis_filename, args.output)
			if args.skip_preprocessing:
				for string in ["gene_predictions_{}_{}.txt","mapped_feature_weights_{}_{}.txt","PSS_{}_{}.txt","GSS_{}_{}.txt"]:
					if os.path.exists(os.path.join(args.output, hypothesis_filename.replace("hypothesis.txt", string.format(lambda_pair1[0], lambda_pair1[1])))):
						os.remove(os.path.join(args.output, hypothesis_filename.replace("hypothesis.txt", string.format(lambda_pair1[0], lambda_pair1[1]))))
				if os.path.exists(hypothesis_filename.replace(".txt","_out_feature_weights_{}_{}.xml".format(lambda_pair1[0], lambda_pair1[1]))) and not args.skip_processing and not args.preserve_xml:
					os.remove(hypothesis_filename.replace(".txt","_out_feature_weights_{}_{}.xml".format(lambda_pair1[0], lambda_pair1[1])))
			if not args.preserve_xml:
				try:
					shutil.move(hypothesis_filename.replace(".txt","_out_feature_weights_{}_{}.xml".format(lambda_pair1[0], lambda_pair1[1])), args.output)
				except:
					pass
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "gene_predictions_{}_{}.txt".format(lambda_pair1[0], lambda_pair1[1])), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "hypothesis_out_feature_weights_{}_{}.txt".format(lambda_pair1[0], lambda_pair1[1])), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "PSS_{}_{}.txt".format(lambda_pair1[0], lambda_pair1[1])), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "GSS_{}_{}.txt".format(lambda_pair1[0], lambda_pair1[1])), args.output)
		#shutil.move(hypothesis_filename.replace("hypothesis.txt", "xval_groups.txt"), args.output)
	# if args.auto_name_nodes:
	#	shutil.move("auto_named_{}".format(os.path.basename(args.tree)), args.output)
	with open(os.path.join(args.output, "HSS.txt"), 'w') as file:
		file.write("{}\t{}\n".format("Hypothesis", "HSS"))
		for hypothesis_filename in hypothesis_file_list:
			file.write("{}\t{}\n".format(hypothesis_filename.replace("_hypothesis.txt", ""), HSS[hypothesis_filename]))
			if args.slep_sample_balance or args.smart_sampling:
				shutil.move(hypothesis_filename.replace("hypothesis.txt", "slep_opts.txt"), args.output)
				shutil.move(hypothesis_filename.replace("hypothesis.txt", "sweights.txt"), args.output)
			if not args.grid_summary_only:
				for lambda_pair1 in lambda_list:
					predictions_fname = hypothesis_filename.replace("hypothesis.txt", "gene_predictions_{}_{}.txt".format(lambda_pair1[0], lambda_pair1[1]))
					if predictions_fname not in missing_predictions:
						try:
							gcv_files.append(gcv.main(os.path.join(args.output,predictions_fname), species_limit=args.species_display_limit, gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff, m_grid=args.m_grid))
						except:
							print("Could not produce graphic for lambda {} with hypothesis {}.".format(lambda_pair1, hypothesis_filename.replace("_hypothesis.txt", "")))
	for file in gcv_files:
		if os.path.dirname(file)!=os.path.normpath(args.output):
			shutil.move(file, args.output)
	return hypothesis_file_list, missing_predictions


def main(args):
	all_files = {"HSS": "HSS.txt"}
	try:
		hypothesis_file_list, slep_opts_file_list = pf.generate_hypothesis_set(args)
		print(hypothesis_file_list)
		print(slep_opts_file_list)
		all_files["hypothesis"] = [os.path.basename(file) for file in hypothesis_file_list]
		all_files["slep_opts"] = [fname.replace("hypothesis.txt", "slep_opts.txt") for fname in all_files["hypothesis"]]
		all_files["sweights"] = [fname.replace("slep_opts", "sweights") for fname in all_files["slep_opts"]]
		all_files["xval"] = [fname.replace("slep_opts", "xval_groups") for fname in all_files["slep_opts"]]
		HSS = {}
		missing_seqs = set()
		merged_rep_predictions_files = {hypothesis_filename:[] for hypothesis_filename in hypothesis_file_list}
		gcv_files = []
		if args.ensemble_parts is not None and args.ensemble_parts >= 1:
			tempdir_list = []
			merged_parts_prediction_files = {hypothesis_filename:[] for hypothesis_filename in hypothesis_file_list}
			for i in range(0, args.ensemble_coverage):
				partitioned_aln_lists = pf.split_gene_list(args.aln_list, args.ensemble_parts)
				j = 0
				gene_prediction_files = {hypothesis_filename:[] for hypothesis_filename in hypothesis_file_list}
				for part_aln_list in partitioned_aln_lists:
					tempdir = "{}_rep{}_part{}".format(args.output, i+1, j+1)
					tempdir_list.append(tempdir)
					features_filename_list, groups_filename_list, response_filename_list, gene_list, field_filename_list, group_list = pf.generate_input_matrices(part_aln_list, hypothesis_file_list, args)
					with open(os.path.join(args.output, "missing_seqs_" + args.output + ".txt"), "r") as file:
						for line in file:
							data = line.strip().split("\t")
							missing_seqs.add((data[1], os.path.splitext(os.path.basename(data[0]))[0]))
					weights_file_list = pf.run_esl(features_filename_list, groups_filename_list, response_filename_list, field_filename_list, args, slep_opts_file_list)
					pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs, group_list)
					for hypothesis_filename in hypothesis_file_list:
						if i+1 == args.ensemble_coverage and j+1 == args.ensemble_parts and not args.sparsify:
							shutil.move(hypothesis_filename, args.output)
						else:
							shutil.copy(hypothesis_filename, args.output)
						shutil.move(hypothesis_filename.replace(".txt","_out_feature_weights.xml"), args.output)
						shutil.move(hypothesis_filename.replace("hypothesis.txt","gene_predictions.txt"), args.output)
						gene_prediction_files[hypothesis_filename].append(os.path.join(tempdir, hypothesis_filename.replace("hypothesis.txt","gene_predictions.txt")))
						shutil.move(hypothesis_filename.replace("hypothesis.txt", "mapped_feature_weights.txt"), args.output)
						shutil.move(hypothesis_filename.replace("hypothesis.txt", "PSS.txt"), args.output)
						shutil.move(hypothesis_filename.replace("hypothesis.txt", "GSS.txt"), args.output)
						shutil.move(hypothesis_filename.replace("hypothesis.txt", "xval_groups.txt"), args.output)
					shutil.move(args.output, tempdir)
					j += 1
				if not args.sparsify:
					for hypothesis_filename in hypothesis_file_list:
						merged_parts_prediction_files[hypothesis_filename].append(gcv.merge_predictions(gene_prediction_files[hypothesis_filename],hypothesis_filename.replace("hypothesis.txt","merged_gene_predictions_rep{}.txt".format(i))))
			if not args.sparsify:
				for hypothesis_filename in hypothesis_file_list:
					merged_rep_predictions_files[hypothesis_filename] = gcv.merge_predictions(merged_parts_prediction_files[hypothesis_filename],hypothesis_filename.replace("hypothesis.txt","merged_gene_predictions_final.txt"))
					gcv_files.extend(merged_parts_prediction_files[hypothesis_filename])
					gcv_files.append(merged_rep_predictions_files[hypothesis_filename])
					gcv_files.append(gcv.main(merged_rep_predictions_files[hypothesis_filename], species_limit=args.species_display_limit, gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff, m_grid=args.m_grid))
			os.mkdir(args.output)
			for tempdir in tempdir_list:
				shutil.move(tempdir, args.output)
			for file in gcv_files:
				shutil.move(file, args.output)
			if args.auto_name_nodes:
				shutil.move("auto_named_{}".format(os.path.basename(args.tree)), args.output)
			with open(os.path.join(args.output, "HSS.txt"), 'w') as file:
				file.write("{}\t{}\n".format("Hypothesis", "HSS"))
				for hypothesis_filename in hypothesis_file_list:
					file.write("{}\t{}\n".format(hypothesis_filename.replace("_hypothesis.txt", ""), HSS[hypothesis_filename]))
					if args.slep_sample_balance:
						shutil.move(hypothesis_filename.replace("hypothesis.txt", "slep_opts.txt"), args.output)
						shutil.move(hypothesis_filename.replace("hypothesis.txt", "sweights.txt"), args.output)
			result_files_list = pf.find_result_files(args, hypothesis_file_list)
			weights = pf.parse_result_files(args, result_files_list)
			return pf.analyze_ensemble_weights(args, weights)
		else:
			features_filename_list, groups_filename_list, response_filename_list, gene_list, field_filename_list, group_list = pf.generate_input_matrices(args.aln_list, hypothesis_file_list, args)
			for key in ["features_filename_list", "groups_filename_list", "response_filename_list", "gene_list", "field_filename_list", "group_list"]:
				if "filename_list" in key:
					all_files[key] = [os.path.basename(file) for file in eval(key)]
			with open(os.path.join(args.output, "missing_seqs_" + args.output + ".txt"), "r") as file:
				for line in file:
					data = line.strip().split("\t")
					missing_seqs.add((data[1], os.path.splitext(os.path.basename(data[0]))[0]))
			hypothesis_list = list(set(["_".join(x.replace("_" + args.output, "").split("_")[0:-1]) for x in hypothesis_file_list]))
			all_files["gene_predictions"] = []
			all_files["feature_weights"] = []
			all_files["GSS"] = []
			all_files["PSS"] = []
			all_files["SPS_SPP"] = []
			all_files["GSS_median"] = []
			all_files["PSS_median"] = []
			all_files["SPS_SPP_median"] = []
			all_files["GCS_median"] = []
			all_files["slep_opts"] = []
			all_files["sweights"] = []
			all_files["xval"] = []
			all_files["pos_stats"] = []
			all_files["response"] = []
			all_files["features"] = []
			all_files["feature_mapping"] = []
			all_files["field"] = []
			all_files["groups"] = []
			all_files["prediction_pngs"] = []
			all_files["roc_pngs"] = []
			all_files["missing_seqs"] = ["missing_seqs_{}.txt".format(args.output)]
			all_files["lambda_list"] = ["{}_lambda_list.txt".format(args.output)]
			for hypothesis in hypothesis_list:
				all_files["GSS_median"] = all_files["GSS_median"] + ["{}_GSS_median.txt".format(hypothesis)]
				all_files["PSS_median"] = all_files["PSS_median"] + ["{}_PSS_median.txt".format(hypothesis)]
				all_files["SPS_SPP_median"] = all_files["SPS_SPP_median"] + ["{}_SPS_SPP_median.txt".format(hypothesis)]
				all_files["GCS_median"] = all_files["GCS_median"] + ["{}_GCS_median.txt".format(hypothesis)]
				all_files["slep_opts"] = all_files["slep_opts"] + ["{}_{}_slep_opts.txt".format(hypothesis, args.output)]
				all_files["sweights"] = all_files["sweights"] + ["{}_{}_sweights.txt".format(hypothesis, args.output)]
				all_files["xval"] = all_files["xval"] + ["{}_{}_xval_groups.txt".format(hypothesis, args.output)]
				all_files["pos_stats"] = all_files["pos_stats"] + ["pos_stats_{}_{}_hypothesis.txt".format(hypothesis, args.output)]
				all_files["response"] = all_files["response"] + ["response_{}_{}_hypothesis.txt".format(hypothesis, args.output)]
				all_files["features"] = all_files["features"] + ["feature_{}_{}_hypothesis.txt".format(hypothesis, args.output)]
				all_files["feature_mapping"] = all_files["feature_mapping"] + ["feature_mapping_{}_{}_hypothesis.txt".format(hypothesis, args.output)]
				all_files["field"] = all_files["field"] + ["field_{}_{}_hypothesis.txt".format(hypothesis, args.output)]
				all_files["groups"] = all_files["groups"] + ["group_indices_{}_{}_hypothesis.txt".format(hypothesis, args.output)]
				all_files["feature_weights"] = all_files["feature_weights"] + ["{}_{}_hypothesis_out_feature_weights.txt".format(hypothesis, args.output)]
				all_files["feature_weights"] = all_files["feature_weights"] + ["{}_{}_hypothesis_out_feature_weights.xml".format(hypothesis, args.output)]
				all_files["gene_predictions"] = all_files["gene_predictions"] + ["{}_{}_gene_predictions.txt".format(hypothesis, args.output)]
				all_files["gene_predictions"] = all_files["gene_predictions"] + ["{}_{}_gene_predictions_xval.txt".format(hypothesis, args.output)]
				all_files["prediction_pngs"] = all_files["prediction_pngs"] + ["{}_{}_gene_predictions.png".format(hypothesis, args.output)]
				all_files["prediction_pngs"] = all_files["prediction_pngs"] + ["{}_{}_gene_predictions_xval.png".format(hypothesis, args.output)]
				all_files["roc_pngs"] = all_files["roc_pngs"] + ["{}_{}_gene_predictions_ROC.png".format(hypothesis, args.output)]
				all_files["SPS_SPP"] = all_files["SPS_SPP"] + ["{}_{}_SPS_SPP.txt".format(hypothesis, args.output)]
				all_files["SPS_SPP"] = all_files["SPS_SPP"] + ["{}_{}_SPS_SPP_xval.txt".format(hypothesis, args.output)]
				all_files["GSS"] = all_files["GSS"] + ["{}_{}_GSS.txt".format(hypothesis, args.output)]
				all_files["PSS"] = all_files["PSS"] + ["{}_{}_PSS.txt".format(hypothesis, args.output)]
			weights_file_list = pf.run_esl(features_filename_list, groups_filename_list, response_filename_list, field_filename_list, args, slep_opts_file_list)
			if args.xval > 1:
				for hypothesis_filename in hypothesis_file_list:
					shutil.copy(os.path.join(args.output, hypothesis_filename.replace("hypothesis.txt", "xval_groups.txt")), hypothesis_filename.replace("hypothesis.txt", "xval_groups.txt"))
	#			pf.process_xval_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, args.xval, missing_seqs, group_list)
				pf.process_sparse_xval_weights(weights_file_list, hypothesis_file_list, groups_filename_list, args.aln_list, gene_list, args.xval, missing_seqs, group_list)
				for hypothesis_filename in hypothesis_file_list:
					os.remove(hypothesis_filename.replace("hypothesis.txt", "xval_groups.txt"))
			# pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs, group_list)
			pf.process_sparse_weights(weights_file_list, hypothesis_file_list, groups_filename_list, HSS, missing_seqs, args)
			for hypothesis_filename in hypothesis_file_list:
				shutil.move(hypothesis_filename, args.output)
				try:
					shutil.move(hypothesis_filename.replace(".txt","_out_feature_weights.xml"), args.output)
				except:
					pass
				shutil.move(hypothesis_filename.replace("hypothesis.txt", "gene_predictions.txt"), args.output)
				if args.xval > 1:
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "gene_predictions_xval.txt"), args.output)
				# shutil.move(hypothesis_filename.replace("hypothesis.txt", "mapped_feature_weights.txt"), args.output)
				# shutil.move(hypothesis_filename.replace("hypothesis.txt", "PSS.txt"), args.output)
				shutil.move(hypothesis_filename.replace("hypothesis.txt", "GSS.txt"), args.output)
				#shutil.move(hypothesis_filename.replace("hypothesis.txt", "xval_groups.txt"), args.output)
			if args.auto_name_nodes:
				shutil.move("auto_named_{}".format(os.path.basename(args.tree)), args.output)
			with open(os.path.join(args.output, "HSS.txt"), 'w') as file:
				file.write("{}\t{}\n".format("Hypothesis", "HSS"))
				for hypothesis_filename in hypothesis_file_list:
					file.write("{}\t{}\n".format(hypothesis_filename.replace("_hypothesis.txt", ""), HSS[hypothesis_filename]))
					if args.slep_sample_balance or args.smart_sampling:
						shutil.move(hypothesis_filename.replace("hypothesis.txt", "slep_opts.txt"), args.output)
						shutil.move(hypothesis_filename.replace("hypothesis.txt", "sweights.txt"), args.output)
					gcv_files.append(gcv.main(os.path.join(args.output,hypothesis_filename.replace("hypothesis.txt", "gene_predictions.txt")), species_limit=args.species_display_limit, gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff, m_grid=args.m_grid))
					if args.xval > 1:
						print("{}-Fold Cross Validation Accuracy:".format(args.xval))
						gcv_files.append(gcv.main(os.path.join(args.output, hypothesis_filename.replace("hypothesis.txt", "gene_predictions_xval.txt")), species_limit=args.species_display_limit, gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff, m_grid=args.m_grid))
			for file in gcv_files:
				if os.path.dirname(file)!=os.path.normpath(args.output):
					shutil.move(file, args.output)
	except Exception as e:
		raise Exception(e)
	finally:
		message_triggered = False
		for key in all_files.keys():
			for file in all_files[key]:
				try:
					shutil.move(file, os.path.join(args.output, file))
					if not message_triggered:
						print("Exception encountered, attempting to clean up files...")
						message_triggered = True
					print("Moved {} to {}.".format(file, args.output))
				except:
					pass
		if single_lambda_pair:
			for summary_type in [output_type for output_type in all_files.keys() if "median" in output_type]:
				for fname in all_files[summary_type]:
					try:
						os.remove(os.path.join(args.output, fname))
					except:
						pass
		for output_type in [mode for mode in ["B", "P", "G", "H"] if mode not in args.stats_out]:
			for fname in all_files.get("{}SS".format(output_type), []):
				try:
					os.remove(os.path.join(args.output, fname))
				except:
					pass
		for output_type in [mode for mode in ["B", "P", "G", "H"] if mode not in args.stats_out]:
			for fname in all_files.get("{}SS_median".format(output_type), []):
				try:
					os.remove(os.path.join(args.output, fname))
				except:
					pass
		if "S" not in args.stats_out:
			for fname in list(all_files.get("SPS_SPP", [])) + list(all_files.get("SPS_SPP_median", [])):
				try:
					os.remove(os.path.join(args.output, fname))
				except:
					pass
		for output_type in ["SPS_SPP", "PSS", "GSS", "hypothesis"]:
			for fname in all_files.get(output_type, []):
				new_fname = "{}_{}".format(output_type, fname.replace("_{}".format(output_type), "").replace("_{}".format(args.output), ""))
				try:
					shutil.move(os.path.join(args.output, fname), os.path.join(args.output, new_fname))
				except:
					pass
		for output_type in ["SPS_SPP_median", "PSS_median", "GSS_median"]:
			for fname in all_files.get(output_type, []):
				new_fname = "{}_{}".format(output_type.replace("_median", ""), fname.replace(output_type, "summary"))
				try:
					shutil.move(os.path.join(args.output, fname), os.path.join(args.output, new_fname))
				except:
					pass
		for fname in all_files.get("feature_weights", []):
			new_fname = "MyESL_model_{}".format(fname.replace("_{}_hypothesis_out_feature_weights".format(args.output), ""))
			try:
				shutil.move(os.path.join(args.output, fname), os.path.join(args.output, new_fname))
			except:
				pass
		for fname in all_files.get("gene_predictions", []):
			new_fname = "GSC_{}".format(fname.replace("_{}_gene_predictions".format(args.output), ""))
			try:
				shutil.move(os.path.join(args.output, fname), os.path.join(args.output, new_fname))
			except:
				pass
		for fname in all_files.get("prediction_pngs", []):
			new_fname = "GSC_{}".format(fname.replace("_{}_gene_predictions".format(args.output), ""))
			try:
				shutil.move(os.path.join(args.output, fname), os.path.join(args.output, new_fname))
			except:
				pass
		for fname in all_files.get("GCS_median", []):
			new_fname = "M-Grid_{}".format(fname.replace("_GCS_median", ""))
			try:
				shutil.move(os.path.join(args.output, fname), os.path.join(args.output, new_fname))
				shutil.move(os.path.join(args.output, fname.replace(".txt", ".png")), os.path.join(args.output, new_fname.replace(".txt", ".png")))
			except:
				pass
		if not args.preserve_inputs:
			for output_type in ["features", "feature_mapping", "field", "pos_stats", "response", "xval", "missing_seqs", "groups", "slep_opts", "sweights"]:
				for fname in all_files.get(output_type, []):
					try:
						os.remove(os.path.join(args.output, fname))
					except:
						pass
	return hypothesis_file_list


def lambda_val2label(lambda_val):
	# lambda_str = "{:g}".format(lambda_val)
	# if lambda_str[-4:] == "0001" and lambda_val > 0.0001:
	# 	lambda_str = lambda_str[0:-4]
	# if lambda_str[0:2] == "0.":
	# 	lambda_str = lambda_str[2:]
	# return lambda_str.rstrip('0')
	if "{:g}".format(lambda_val)[0:2] == "0.":
		return "{:g}".format(lambda_val)[2:]
	else:
		return "{:g}".format(lambda_val)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Phylogenetic hypothesis tester.")
	parser.add_argument("aln_list", help="List of alignment files to extract features from.", type=str)
	#parser.add_argument("--response", help="File containing list of named node/response value pairs.", type=str, default=None)
	parser.add_argument("--classes", help="File containing list of named node/response value pairs.", type=str, default=None)
	parser.add_argument("--tree", help="Input phylogeny to perform testing on.", type=str, default=None)
	#parser.add_argument("--nodelist", help="File containing list of named internal nodes to test. If no file is specified, each named internal node of input phylogeny will be tested.", type=str, default=None)
	parser.add_argument("--clade_list", help="File containing list of named internal nodes to test. If no file is specified, each named internal node of input phylogeny will be tested.", type=str, default=None)
	#parser.add_argument("--auto_name_nodes", help="Assign automatically generated names to unnamed internal nodes, causing them to be tested.", action='store_true', default=False)
	parser.add_argument("--gen_clade_list", help="<lower_bound,upper_bound> Assign automatically generated names to unnamed internal nodes containing between lower_bound and upper_bound number of species, causing them to be tested.", type=str, default=None)
	parser.add_argument("--auto_name_length", help="Number of characters to take from sequence IDs to generate internal node labels.", type=int, default=5)
	parser.add_argument("--cladesize_cutoff_lower", help="Internal nodes with fewer than cladesize_cutoff_lower terminal descendants will not be tested.", type=int, default=0)
	parser.add_argument("--cladesize_cutoff_upper", help="Internal nodes with greater than cladesize_cutoff_upper terminal descendants will not be tested.", type=int,default=None)
	# parser.add_argument("--smart_sampling", help="For each selected positive sample, select a balanced, phylogenetically informed negative sample.", type=int, default=None)
	# parser.add_argument("--slep_sample_balance", help="Automatically uses the SLEP option to balance sample weights, only available for non-overlapping, logistic sglasso.",
	# 					action='store_true', default=False)
	# parser.add_argument("--upsample_balance", help="Balance positive and negative response sets by upsampling the underpopulated set.", action='store_true', default=False)
	# parser.add_argument("--downsample_balance", help="Balance positive and negative response sets by downsampling the overpopulated set.", action='store_true', default=False)
	parser.add_argument("--class_bal", help="Sample balancing type:['weighted', 'up', 'down', 'phylo'(, 'phylo_1', 'phylo_2')]", type=str, default=None)
	parser.add_argument("--skip_preprocessing", help="Assume preprocessing files have already been generated.", action='store_true', default=False)
	parser.add_argument("--skip_processing", help="Assume results XML files have already been generated.", action='store_true', default=False)
	parser.add_argument("--preserve_xml", help="Leave results XML files in place for reprocessing.", action='store_true', default=False)
	parser.add_argument("--preserve_inputs", help="Leave input files in place for reprocessing.", action='store_true', default=False)
	parser.add_argument("-z", "--lambda1", help="Feature sparsity parameter.", type=float, default=0.1)
	parser.add_argument("-y", "--lambda2", help="Group sparsity parameter.", type=float, default=0.1)
	#parser.add_argument("--grid_z", help="Grid search sparsity parameter interval specified as 'min,max,step_size'", type=str, default=None)
	#parser.add_argument("--grid_y", help="Grid search group sparsity parameter interval specified as 'min,max,step_size'", type=str, default=None)
	parser.add_argument("--lambda1_range", help="Grid search sparsity parameter interval specified as 'min,max,step_size'", type=str, default=None)
	parser.add_argument("--lambda2_range", help="Grid search group sparsity parameter interval specified as 'min,max,step_size'", type=str, default=None)
	parser.add_argument("--grid_rmse_cutoff", help="RMSE cutoff when selecting models to aggregate.", type=float, default=100.0)
	parser.add_argument("--grid_acc_cutoff", help="Accuracy cutoff when selecting models to aggregate.", type=float, default=0.0)
	parser.add_argument("--grid_threads", help="Number of threads to use when aggregating grid search results.", type=int, default=None)
	parser.add_argument("--grid_summary_only", help="Skip generating graphics for individual runs/models.", action='store_true', default=False)
	#parser.add_argument("--grid_gene_threshold", help="Stop increasing L2 when model selects this many genes or fewer.", type=int, default=None)
	parser.add_argument("--min_groups", help="Stop increasing L2 when model selects this many genes or fewer.", type=int, default=None)
	parser.add_argument("--no_group_penalty", help="Perform mono-level optimization, ignoring group level sparsity penalties.",
						action='store_true', default=False)
	parser.add_argument("-o", "--output", help="Output directory.", type=str, default="output")
	parser.add_argument("--fuzz_indels", help="Assign indels in variable sites a one-hot value of 0.5 instead of 0.", action='store_true', default=False)
	parser.add_argument("--ensemble_parts", help="Build gene-wise ensemble models, splitting the set of genes into N partitions for each run.", type=int, default=None)
	parser.add_argument("--ensemble_coverage", help="Number of ensemble models to build. Each gene will be included in this many individual models.", type=int, default=5)
	parser.add_argument("--sparsify", help="Iteratively increase sparsity until selected set of genes fits in one partition.", action='store_true', default=False)
	parser.add_argument("--method", help="SGLasso type to use. Options are \"logistic\", \"leastr\", or \"ol_leastr\". Defaults to \"leastr\".", type=str, default="logistic")
	parser.add_argument("--slep_opts", help="File of tab-separated name-value pairs (one per line) to specify SLEP options.", type=str, default=None)
	#parser.add_argument("--gene_penalties", help="File of penalty values (same order as aln_list) to specify penalty score for each gene.", type=str, default=None)
	parser.add_argument("--group_wt", help="File of penalty values (same order as aln_list) to specify penalty score for each gene.", type=str, default=None)
	#parser.add_argument("--xval", help="Number of partitions to use when performing cross-validation.", type=int, default=1)
	parser.add_argument("--kfold", help="Number of partitions to use when performing cross-validation.", type=int, default=1)
	parser.add_argument("--gene_display_limit", help="Limits the number of genes displayed in the generated graph images.", type=int, default=100)
	parser.add_argument("--gene_display_cutoff", help="Limits genes displayed in the generated graph images to those with sum-of-squares greater than cutoff value.", type=int, default=0.0)
	parser.add_argument("--species_display_limit", help="Limits the number of species displayed in the generated graph images.", type=int, default=100)
	#parser.add_argument("--m_grid", help="Generate m-grid graphical output.", action='store_true', default=False)
	parser.add_argument("--m_grid", help="Generate m-grid graphical output with display limited to <rows>,<cols>, e.g. \"--m_grid 10,20\" .", type=str, default=None)
	parser.add_argument("--bit_ct", help="Ignore mutations observed fewer than N times.", type=int, default=1)
	parser.add_argument("--stats_out", help="<str[BPGHS]*> Various single-character flags for output produced, consult README for more info.", type=str, default="")
	parser.add_argument("--data_type", help="<Options are \"nucleotide\", \"protein\", \"molecular\", \"universal\". Consult documentation for detailed info.", type=str, default="universal")
	parser.add_argument("--DrPhylo", help="Run Dr Phylo type analysis.", action='store_true', default=False)
	args = parser.parse_args()
	if args.classes is None and args.tree is None:
		raise Exception("Must invoke one of --tree or --classes option.")
	if args.classes is not None and args.tree is not None:
		raise Exception("Cannot use --tree and --classes options simultaneously.")
	if args.tree is None and args.gen_clade_list:
		raise Exception("Cannot use --gen_clade_list option and --classes options simultaneously.")
	if args.DrPhylo:
		if args.class_bal is None:
			if args.classes is not None:
				args.class_bal = "weighted"
			if args.tree is not None:
				args.class_bal = "phylo"
		if args.lambda1_range is None:
			args.lambda1_range = "0.1,1.0,0.1"
		if args.lambda2_range is None:
			args.lambda2_range = "0.1,1.0,0.1"
		if args.min_groups is None:
			args.min_groups = 1
		if args.stats_out == "":
			args.stats_out = "BPGHS"
		if args.m_grid is None:
			args.m_grid = "20,30"
	if args.lambda1_range is not None and args.lambda2_range is not None and args.kfold > 1:
		raise Exception("Cannot use --kfold option while running in grid search mode.")
	single_lambda_pair = False
	#Convert single run to 1x1 grid run
	if args.lambda1_range is None and args.lambda2_range is None:
		args.lambda1_range = "{},{},{}".format(args.lambda1, args.lambda1 + 0.00001, 0.00002)
		args.lambda2_range = "{},{},{}".format(args.lambda2, args.lambda2 + 0.00001, 0.00002)
		single_lambda_pair = True
	### Convert new parameter names to old ###
	args.slep_sample_balance = False
	args.upsample_balance = False
	args.downsample_balance = False
	args.smart_sampling = None
	args.response = args.classes
	args.xval = args.kfold
	args.grid_gene_threshold = args.min_groups
	if args.class_bal == "weighted":
		args.slep_sample_balance = True
	elif args.class_bal == "up":
		args.upsample_balance = True
	elif args.class_bal == "down":
		args.downsample_balance = True
	elif args.class_bal == "phylo":
		args.smart_sampling = 2
	elif args.class_bal == "phylo_1":
		args.smart_sampling = 1
	elif args.class_bal == "phylo_2":
		args.smart_sampling = 3
	args.grid_z = args.lambda1_range
	args.grid_y = args.lambda2_range
	args.auto_name_nodes = False
	if args.gen_clade_list is not None:
		args.auto_name_nodes = True
		cutoffs = args.gen_clade_list.strip().split(",")
		args.cladesize_cutoff_lower = int(cutoffs[0])
		args.cladesize_cutoff_upper = int(cutoffs[1])
	args.nodelist = args.clade_list
	args.gene_penalties = args.group_wt
	if args.m_grid is None:
		args.m_grid = False
	else:
		m_grid_dims = args.m_grid.strip().split(",")
		args.species_display_limit = int(m_grid_dims[0])
		args.gene_display_limit = int(m_grid_dims[1])
		args.m_grid = True
	if args.no_group_penalty:
		args.lambda2 = 0.0000000001
	if (args.grid_z is None and args.grid_y is not None) or (args.grid_y is None and args.grid_z is not None):
		raise Exception("Only one grid search parameter specified, --grid_z and --grid_y must be specified together.")
	elif args.grid_z is not None and args.grid_y is not None and args.xval <= 1:
		all_files = {"HSS": "HSS.txt"}
		try:
			args_original = copy.deepcopy(args)
			output_folder = args.output
			if not os.path.exists(output_folder):
				os.mkdir(output_folder)
			hypothesis_file_list = []
			# for x in args.grid_z.strip().split(','):
			# 	print(x)
			# 	print(float(x))
			[z_min, z_max, z_interval] = [float(x) for x in args.grid_z.strip().split(',')]
			#z_interval = (z_max - z_min) / z_steps
			z_steps = (z_max - z_min) / z_interval
			z_list = [z_min + (x * z_interval) for x in range(0, int(z_steps) + 1) if z_min + (x * z_interval) < 1.0]
			[y_min, y_max, y_interval] = [float(x) for x in args.grid_y.split(',')]
			#y_interval = (y_max - y_min) / y_steps
			y_steps = (y_max - y_min) / y_interval
			y_list = [y_min + (x * y_interval) for x in range(0, int(y_steps) + 1) if y_min + (x * y_interval) < 1.0]
			[z_count, y_count] = [0, 0]
			start = datetime.now()
			input_files = {}
			input_files["hypothesis_file_list"], input_files["slep_opts_file_list"] = pf.generate_hypothesis_set(args)
			print(input_files["slep_opts_file_list"])
			all_files["hypothesis"] = [os.path.basename(file) for file in input_files["hypothesis_file_list"]]
			all_files["slep_opts"] = [os.path.basename(file) for file in input_files["slep_opts_file_list"]]
			all_files["sweights"] = [fname.replace("slep_opts", "sweights") for fname in all_files["slep_opts"]]
			all_files["xval"] = [fname.replace("slep_opts", "xval_groups") for fname in all_files["slep_opts"]]
#			raise Exception("Testing")
			input_files["features_filename_list"], input_files["groups_filename_list"], input_files["response_filename_list"], input_files["gene_list"], \
				input_files["field_filename_list"], input_files["group_list"] = pf.generate_input_matrices(args.aln_list, input_files["hypothesis_file_list"], args)
			for key in input_files.keys():
				if "filename_list" in key:
					all_files[key] = [os.path.basename(file) for file in input_files[key]]
			input_files["lambda_list"] = os.path.join(args.output, "{}_lambda_list.txt".format(args.output))
			lambda_list = []
			timers["pre_processing"] = datetime.now() - start
			print("Time elapsed while processing alignments to features: {}".format(timers["pre_processing"]))
			with open(os.path.join(args.output, "{}_lambda_list.txt".format(args.output)), 'w') as file:
				for lambda1 in z_list:
					for lambda2 in y_list:
						lambda_list.append((lambda1, lambda2))
						file.write("{:g}\t{:g}\n".format(lambda1, lambda2))
			# os.mkdir(args_original.output)
			hypothesis_list = list(set(["_".join(x.replace("_" + args_original.output, "").split("_")[0:-1]) for x in input_files["hypothesis_file_list"]]))
			all_files["gene_predictions"] = []
			all_files["feature_weights"] = []
			all_files["GSS"] = []
			all_files["PSS"] = []
			all_files["SPS_SPP"] = []
			all_files["GSS_median"] = []
			all_files["PSS_median"] = []
			all_files["SPS_SPP_median"] = []
			all_files["GCS_median"] = []
			all_files["slep_opts"] = []
			all_files["sweights"] = []
			all_files["xval"] = []
			all_files["pos_stats"] = []
			all_files["response"] = []
			all_files["features"] = []
			all_files["feature_mapping"] = []
			all_files["field"] = []
			all_files["groups"] = []
			all_files["prediction_pngs"] = []
			all_files["roc_pngs"] = []
			all_files["missing_seqs"] = ["missing_seqs_{}.txt".format(args_original.output)]
			all_files["lambda_list"] = ["{}_lambda_list.txt".format(args_original.output)]
			for hypothesis in hypothesis_list:
				all_files["GSS_median"] = all_files["GSS_median"] + ["{}_GSS_median.txt".format(hypothesis)]
				all_files["PSS_median"] = all_files["PSS_median"] + ["{}_PSS_median.txt".format(hypothesis)]
				all_files["SPS_SPP_median"] = all_files["SPS_SPP_median"] + ["{}_SPS_SPP_median.txt".format(hypothesis)]
				all_files["GCS_median"] = all_files["GCS_median"] + ["{}_GCS_median.txt".format(hypothesis)]
				all_files["slep_opts"] = all_files["slep_opts"] + ["{}_{}_slep_opts.txt".format(hypothesis, args.output)]
				all_files["sweights"] = all_files["sweights"] + ["{}_{}_sweights.txt".format(hypothesis, args.output)]
				all_files["xval"] = all_files["xval"] + ["{}_{}_xval_groups.txt".format(hypothesis, args.output)]
				all_files["pos_stats"] = all_files["pos_stats"] + ["pos_stats_{}_{}_hypothesis.txt".format(hypothesis, args_original.output)]
				all_files["response"] = all_files["response"] + ["response_{}_{}_hypothesis.txt".format(hypothesis, args_original.output)]
				all_files["features"] = all_files["features"] + ["feature_{}_{}_hypothesis.txt".format(hypothesis, args_original.output)]
				all_files["feature_mapping"] = all_files["feature_mapping"] + ["feature_mapping_{}_{}_hypothesis.txt".format(hypothesis, args_original.output)]
				all_files["field"] = all_files["field"] + ["field_{}_{}_hypothesis.txt".format(hypothesis, args_original.output)]
				all_files["groups"] = all_files["groups"] + ["group_indices_{}_{}_hypothesis.txt".format(hypothesis, args_original.output)]
				for lambda_pair in lambda_list:
					all_files["feature_weights"] = all_files["feature_weights"] + ["{}_{}_hypothesis_out_feature_weights_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))]
					all_files["feature_weights"] = all_files["feature_weights"] + ["{}_{}_hypothesis_out_feature_weights_{}_{}.xml".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))]
					all_files["gene_predictions"] = all_files["gene_predictions"] + ["{}_{}_gene_predictions_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))]
					all_files["prediction_pngs"] = all_files["prediction_pngs"] + ["{}_{}_gene_predictions_{}_{}.png".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))]
					all_files["roc_pngs"] = all_files["roc_pngs"] + ["{}_{}_gene_predictions_{}_{}_ROC.png".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))]
					all_files["SPS_SPP"] = all_files["SPS_SPP"] + ["{}_{}_SPS_SPP_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))]
					all_files["GSS"] = all_files["GSS"] + ["{}_{}_GSS_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))]
					all_files["PSS"] = all_files["PSS"] + ["{}_{}_PSS_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))]
			hypothesis_file_list, missing_predictions = grid_search(args, args.output, input_files)
			pss_vals = {}
			gss_vals = {}
			# gene_predictions = []
			# print(missing_predictions)
			for hypothesis in hypothesis_list:
				gene_predictions = []
				for lambda_pair in lambda_list:
					# print("{}_{}_gene_predictions_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1])))
					if "{}_{}_gene_predictions_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1])) in missing_predictions:
						continue
					# Parse gene_predictions files
					rmse = 0
					acc = 1
					# model = pandas.read_csv(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind), "{}_{}_{}_{}_gene_predictions.txt".format(hypothesis, args_original.output, z_ind, y_ind)), delim_whitespace=True, index_col=0)
					model = pandas.read_csv(os.path.join(args_original.output, "{}_{}_gene_predictions_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))), delim_whitespace=True, index_col=0)
					rmse = ((model.Response - model.Prediction) ** 2).mean() ** 0.5
					acc = sum([1 for x in zip(model.Response, model.Prediction) if (x[0] == -1 and x[1] < 0) or (x[0] == 1 and x[1] > 0)])/float(model.shape[0])
					# Calculate acc and RMSE
					if rmse <= args_original.grid_rmse_cutoff and acc >= args_original.grid_acc_cutoff:
						gene_predictions.append(model)
					# Parse GSS file into GSS vals
					# with open(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind), "{}_{}_{}_{}_GSS.txt".format(hypothesis, args_original.output, z_ind, y_ind)), 'r') as gss_file:
					with open(os.path.join(args_original.output, "{}_{}_GSS_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))), 'r') as gss_file:
						for line in gss_file:
							data = line.strip().split('\t')
							if data[0] == "Gene":
								continue
							if data[0] in gss_vals:
								gss_vals[data[0]].update({"{}_{}".format(lambda_pair[0], lambda_pair[1]): data[1]})
							else:
								gss_vals[data[0]] = {"{}_{}".format(lambda_pair[0], lambda_pair[1]): data[1]}
					# Parse PSS file into PSS vals
					# with open(os.path.join(args_original.output, args_original.output + "_{}_{}".format(z_ind, y_ind), "{}_{}_{}_{}_PSS.txt".format(hypothesis, args_original.output, z_ind, y_ind)), 'r') as pss_file:
					with open(os.path.join(args_original.output, "{}_{}_PSS_{}_{}.txt".format(hypothesis, args_original.output, lambda_val2label(lambda_pair[0]), lambda_val2label(lambda_pair[1]))), 'r') as pss_file:
						for line in pss_file:
							data = line.strip().split('\t')
							if data[0] == "Position Name" or float(data[1]) == 0.0:
								continue
							if data[0] in pss_vals:
								pss_vals[data[0]].update({"{}_{}".format(lambda_pair[0], lambda_pair[1]): data[1]})
							else:
								pss_vals[data[0]] = {"{}_{}".format(lambda_pair[0], lambda_pair[1]): data[1]}
				# Write aggregated GSS file
				with open(os.path.join(args_original.output, "{}_GSS_median.txt".format(hypothesis)), 'w') as file:
					for gene in gss_vals.keys():
						gss_temp = [float(x) for x in gss_vals[gene].values() if abs(float(x)) > 0]
						if len(gss_temp) == 0:
							gss_temp = [0]
						file.write("{}\t{}\n".format(gene, statistics.median(gss_temp)))
				# Write aggregated PSS file
				with open(os.path.join(args_original.output, "{}_PSS_median.txt".format(hypothesis)), 'w') as file:
					for pos in pss_vals.keys():
						pss_temp = [float(x) for x in pss_vals[pos].values() if abs(float(x)) > 0]
						if len(pss_temp) == 0:
							pss_temp = [0]
						file.write("{}\t{}\n".format(pos, statistics.median(pss_temp)))
				# Write aggregated gene_predictions file
				with open(os.path.join(args_original.output, "{}_GCS_median.txt".format(hypothesis)), 'w') as file:
					gene_list = []
					species_list = set()
					for gene_prediction in gene_predictions:
						species_list.update(gene_prediction.index)
						for column in gene_prediction.columns:
							if column not in gene_list:
								gene_list.append(column)
					file.write("SeqID\tResponse\tPrediction_mean\t{}\n".format("\t".join([gene for gene in gene_list if gene not in ["SeqID", "Prediction", "Intercept", "Response"]])))
					species_list = list(species_list)
					#for species in gene_predictions[0].index:
					for species in species_list:
						file.write("{}\t".format(species))
						for gene in gene_list:
							if gene in ["SeqID", "Intercept"]:
								continue
							GCS_list = [x[gene][species] for x in gene_predictions if gene in x.columns.values.tolist() and species in x[gene].keys() and x[gene][species] != 0]
							if len(GCS_list) == 0:
								GCS_list = [0]
							# if gene == "Response":
							# 	print(GCS_list)
							# 	print(statistics.median(GCS_list))
							if gene == "Prediction":
								#file.write("{}\t".format(statistics.mode(GCS_list)))
								file.write("{}\t".format(statistics.mean(GCS_list)))
							else:
								file.write("{}\t".format(statistics.median(GCS_list)))
						file.write("\n")
				gcv.main(os.path.join(args_original.output, "{}_GCS_median.txt".format(hypothesis)), lead_cols=3, species_limit=args.species_display_limit, gene_limit=args.gene_display_limit, ssq_threshold=args.gene_display_cutoff, m_grid=args.m_grid)
			for file in input_files["hypothesis_file_list"]:
				if args.skip_preprocessing and os.path.exists(os.path.join(args_original.output, file)):
					os.remove(os.path.join(args_original.output, file))
				shutil.move(file, args_original.output)
		except Exception as e:
			raise Exception(e)
		finally:
			message_triggered = False
			for key in all_files.keys():
				for file in all_files[key]:
					try:
						shutil.move(file, os.path.join(args_original.output, file))
						if not message_triggered:
							print("Exception encountered, attempting to clean up files...")
							message_triggered = True
						print("Moved {} to {}.".format(file, args_original.output))
					except:
						pass
			if single_lambda_pair:
				for summary_type in [output_type for output_type in all_files.keys() if "median" in output_type]:
					for fname in all_files[summary_type]:
						try:
							os.remove(os.path.join(args_original.output, fname))
						except:
							pass
			for output_type in [mode for mode in ["B", "P", "G", "H"] if mode not in args.stats_out]:
				for fname in all_files.get("{}SS".format(output_type), []):
					try:
						os.remove(os.path.join(args_original.output, fname))
					except:
						pass
			for output_type in [mode for mode in ["B", "P", "G", "H"] if mode not in args.stats_out]:
				for fname in all_files.get("{}SS_median".format(output_type), []):
					try:
						os.remove(os.path.join(args_original.output, fname))
					except:
						pass
			if "S" not in args.stats_out:
				for fname in list(all_files.get("SPS_SPP", [])) + list(all_files.get("SPS_SPP_median", [])):
					try:
						os.remove(os.path.join(args_original.output, fname))
					except:
						pass
			for output_type in ["SPS_SPP", "PSS", "GSS", "hypothesis"]:
				for fname in all_files.get(output_type, []):
					new_fname = "{}_{}".format(output_type, fname.replace("_{}".format(output_type), "").replace("_{}".format(args_original.output), ""))
					try:
						shutil.move(os.path.join(args_original.output, fname), os.path.join(args_original.output, new_fname))
					except:
						pass
			for output_type in ["SPS_SPP_median", "PSS_median", "GSS_median"]:
				for fname in all_files.get(output_type, []):
					new_fname = "{}_{}".format(output_type.replace("_median", ""), fname.replace(output_type, "summary"))
					try:
						shutil.move(os.path.join(args_original.output, fname), os.path.join(args_original.output, new_fname))
					except:
						pass
			for fname in all_files.get("feature_weights", []):
				new_fname = "MyESL_model_{}".format(fname.replace("_{}_hypothesis_out_feature_weights".format(args_original.output), ""))
				try:
					shutil.move(os.path.join(args_original.output, fname), os.path.join(args_original.output, new_fname))
				except:
					pass
			for fname in all_files.get("gene_predictions", []):
				new_fname = "GSC_{}".format(fname.replace("_{}_gene_predictions".format(args_original.output), ""))
				try:
					shutil.move(os.path.join(args_original.output, fname), os.path.join(args_original.output, new_fname))
				except:
					pass
			for fname in all_files.get("prediction_pngs", []):
				new_fname = "GSC_{}".format(fname.replace("_{}_gene_predictions".format(args_original.output), ""))
				try:
					shutil.move(os.path.join(args_original.output, fname), os.path.join(args_original.output, new_fname))
				except:
					pass
			for fname in all_files.get("GCS_median", []):
				new_fname = "M-Grid_{}".format(fname.replace("_GCS_median", ""))
				try:
					shutil.move(os.path.join(args_original.output, fname), os.path.join(args_original.output, new_fname))
					shutil.move(os.path.join(args_original.output, fname.replace(".txt", ".png")), os.path.join(args_original.output, new_fname.replace(".txt", ".png")))
				except:
					pass
			if not args.preserve_inputs:
				for output_type in ["features", "feature_mapping", "field", "pos_stats", "response", "xval", "missing_seqs", "groups", "slep_opts", "sweights"]:
					for fname in all_files.get(output_type, []):
						try:
							os.remove(os.path.join(args_original.output, fname))
						except:
							pass


	else:
		score_tables = main(args)
	gene_target = None
	if args.sparsify and score_tables is not None:
		args_original = copy.deepcopy(args)
		output_folder = args.output
		aln_list_dir = os.path.dirname(args.aln_list)
		aln_list_basename = os.path.splitext(os.path.basename(args.aln_list))[0]
		aln_file_list = {}
		with open(args.aln_list, 'r') as file:
			for line in file:
				for aln_filename in line.strip().split(","):
					basename = os.path.splitext(os.path.basename(aln_filename.strip()))[0]
					if basename in aln_file_list.keys():
						if aln_file_list[basename] != aln_filename.strip():
							raise Exception("Found multiple alignment files with identical basename {}.".format(basename))
					else:
						aln_file_list[basename] = aln_filename.strip()
		for hypothesis in score_tables.keys():
			args = copy.deepcopy(args_original)
			if gene_target is None:
				gene_target = math.ceil(len(score_tables[hypothesis])/args.ensemble_parts)
			elif gene_target != math.ceil(len(score_tables[hypothesis])/args.ensemble_parts):
				raise Exception("Mismatched gene counts between hypotheses.")
			selected_count = sum([1 for val in score_tables[hypothesis].values() if val[0] > 0])
			counter = 1
			final_round = False
			if selected_count <= gene_target and not final_round:
				final_round = True
			shutil.move("{}_hypothesis.txt".format(hypothesis), "{}.txt".format(hypothesis))
			while selected_count > gene_target or final_round:
				# Generate new aln_list file and point args.aln_list at it
				aln_list_filename = os.path.join(aln_list_dir, "{}_{}_{}.txt".format(aln_list_basename, hypothesis, counter))
				with open(aln_list_filename, 'w') as file:
					with open(args.aln_list, 'r') as original_file:
						stashed_gene_groups = []
						nonzero_genes = [os.path.basename(key) for key in score_tables[hypothesis].keys() if score_tables[hypothesis][key][0] > 0]
						#print(nonzero_genes)
						new_gene_groups = []
						for line in original_file:
							gene_group = line.strip().split(",")
							#print(len(gene_group))
							new_gene_group = sorted([os.path.splitext(os.path.basename(val))[0] for val in gene_group if os.path.splitext(os.path.basename(val))[0] in nonzero_genes])
							print("{}\t::\t{}".format(",".join(gene_group), ",".join(new_gene_group)))
							if tuple(new_gene_group) not in stashed_gene_groups:
								stashed_gene_groups.append(tuple(new_gene_group))
								new_gene_groups.append(new_gene_group)
						#print(new_gene_groups)
						for new_gene_group in new_gene_groups:
							if len(new_gene_group) > 0:
								file.write("{}\n".format(",".join([aln_file_list[val] for val in new_gene_group])))
				args.aln_list = aln_list_filename
				# Also point --response at this hypothesis' response file
				args.response = "{}.txt".format(hypothesis)
				# args.output = os.path.join(output_folder, "{}_sparsify_round{}".format(hypothesis, counter))
				args.output = "{}_sparsify_round{}".format(hypothesis, counter)
				# temp_score_tables = main(args)
				if final_round:
					args.output = "{}_sparsify_final".format(hypothesis)
					args.lambda1 = 10**-6
					args.lambda2 = 10**-6
					args.ensemble_parts = 1
					args.ensemble_coverage = 1
					args.sparsify = False
				score_tables.update(main(args))
				args.sparsify = True
				shutil.move(args.output, output_folder)
				new_selected_count = sum([1 for val in score_tables[hypothesis].values() if val[0] > 0])
				if not final_round and new_selected_count == selected_count:
					# This round didn't increase sparsity, so increase group sparsity parameter for the next round
					print("The number of selected genes was unchanged, increasing group sparsity parameter from {} to {}.".format(args.lambda2, args.lambda2 * 1.2))
					args.lambda2 = args.lambda2 * 1.2
				else:
					selected_count = new_selected_count
				if selected_count <= gene_target and not final_round:
					final_round = True
				else:
					final_round = False
				counter += 1
			if os.path.exists("{}.txt".format(hypothesis)):
				os.remove("{}.txt".format(hypothesis))
			try:
				shutil.move(os.path.join(args_original.output,"{}_sparsify_final".format(hypothesis),"{}_merged_gene_predictions_final.txt".format(hypothesis)), args_original.output)
				shutil.move(os.path.join(args_original.output,"{}_sparsify_final".format(hypothesis),"{}_merged_gene_predictions_final.png".format(hypothesis)), args_original.output)
			except:
				print("Couldn't find summary graphic file at {}.".format(os.path.join(args_original.output,"{}_sparsify_final".format(hypothesis),"{}_merged_gene_predictions_final.png".format(hypothesis))))

	# generate_gene_prediction_table(weights_filename, responses_filename, groups_filename, features_filename, output_filename)
	# if False:
	# 	pf.generate_gene_prediction_table("angiosperm_out_feature_weights.xml", "angiosperm_20spec_pred.txt", "angiosperm_input/group_indices_angiosperm_input.txt", "angiosperm_input/feature_angiosperm_input.txt", "test_output.txt")
	# if False:
	# 	hypothesis_file_list = pf.generate_hypothesis_set("testtree2.nwk", "nodelist.txt")
	# 	features_filename, groups_filename, response_filename_list = pf.generate_input_matrices("sample_files/angiosperm_100_sample_alns.txt", hypothesis_file_list)
	# 	weights_file_list = pf.run_esl(features_filename, groups_filename, response_filename_list)
	# 	pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename, features_filename)
	for key in timers.keys():
		print("Total time elapsed for {}: {}".format(key, timers[key]))

