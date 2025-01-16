import os
import sys
import time
import shutil
import argparse
import math
import copy
import pipeline_funcs as pf






if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Phylogenetic hypothesis tester.")
	parser.add_argument("aln_list", help="List of alignment files to extract features from.", type=str)
	parser.add_argument("--classes", help="File containing list of named node/response value pairs.", type=str, default=None)
	parser.add_argument("--tree", help="Input phylogeny to perform testing on.", type=str, default=None)
	parser.add_argument("--clade_list", help="File containing list of named internal nodes to test. If no file is specified, each named internal node of input phylogeny will be tested.", type=str, default=None)
	parser.add_argument("--gen_clade_list", help="<lower_bound,upper_bound> Assign automatically generated names to unnamed internal nodes containing between lower_bound and upper_bound number of species, causing them to be tested.", type=str, default=None)
	parser.add_argument("--auto_name_length", help="Number of characters to take from sequence IDs to generate internal node labels.", type=int, default=5)
	parser.add_argument("--cladesize_cutoff_lower", help="Internal nodes with fewer than cladesize_cutoff_lower terminal descendants will not be tested.", type=int, default=0)
	parser.add_argument("--cladesize_cutoff_upper", help="Internal nodes with greater than cladesize_cutoff_upper terminal descendants will not be tested.", type=int,default=None)
	parser.add_argument("--class_bal", help="Sample balancing type:['weighted', 'up', 'down', 'phylo'(, 'phylo_1', 'phylo_2')]", type=str, default=None)
	parser.add_argument("--preserve_inputs", help="Leave input files in place for inspection.", action='store_true', default=False)
	parser.add_argument("--disable_mc", help="Disable available memory check.", action='store_true', default=False)
	parser.add_argument("-z", "--lambda1", help="Feature sparsity parameter.", type=float, default=None)
	parser.add_argument("-y", "--lambda2", help="Group sparsity parameter.", type=float, default=None)
	parser.add_argument("--lambda1_grid", help="Grid search sparsity parameter interval specified as 'min,max,step_size'", type=str, default=None)
	parser.add_argument("--lambda2_grid", help="Grid search group sparsity parameter interval specified as 'min,max,step_size'", type=str, default=None)
	parser.add_argument("--grid_rmse_cutoff", help="RMSE cutoff when selecting models to aggregate.", type=float, default=100.0)
	parser.add_argument("--grid_acc_cutoff", help="Accuracy cutoff when selecting models to aggregate.", type=float, default=0.0)
	parser.add_argument("--grid_summary_only", help="Skip generating graphics for individual runs/models and remove individual model files.", action='store_true', default=False)
	parser.add_argument("--min_groups", help="Stop increasing L2 when model selects this many genes or fewer.", type=int, default=-1)
	parser.add_argument("--no_group_penalty", help="Perform mono-level optimization, ignoring group level sparsity penalties.",
						action='store_true', default=False)
	parser.add_argument("-o", "--output", help="Output directory.", type=str, default="output")
	parser.add_argument("--slep_opts", help="File of tab-separated name-value pairs (one per line) to specify SLEP options.", type=str, default=None)
	parser.add_argument("--group_wt", help="File of penalty values (same order as aln_list) to specify penalty score for each gene.", type=str, default=None)
	parser.add_argument("--kfold", help="Number of partitions to use when performing cross-validation.", type=int, default=1)
	parser.add_argument("--kfold_ids", help="File of cross-validation group IDs to reuse.", type=str, default=None)
	parser.add_argument("--gene_display_limit", help="Limits the number of genes displayed in the generated graph images.", type=int, default=100)
	parser.add_argument("--gene_display_cutoff", help="Limits genes displayed in the generated graph images to those with sum-of-squares greater than cutoff value.", type=int, default=0.0)
	parser.add_argument("--species_display_limit", help="Limits the number of species displayed in the generated graph images.", type=int, default=100)
	parser.add_argument("--m_grid", help="Generate m-grid graphical output with display limited to <rows>,<cols>, e.g. \"--m_grid 10,20\" .", type=str, default=None)
	parser.add_argument("--bit_ct", help="Ignore mutations observed fewer than N times.", type=int, default=1)
	parser.add_argument("--auto_bit_ct", help="Automatically set bit_ct to X% of response class with fewer members.", type=float, default=None)
	parser.add_argument("--include_singletons", help="Include singleton sites in analysis.", action='store_true', default=False)
	parser.add_argument("--stats_out", help="<str[PGHS]*> Various single-character flags for output produced, consult README for more info.", type=str, default="")
	parser.add_argument("--data_type", help="<Options are \"nucleotide\", \"protein\", \"molecular\", \"universal\". Consult documentation for detailed info.", type=str, default="universal")
	parser.add_argument("--method", help="SGLasso type to use. Options are \"logistic\" or \"leastr\". Defaults to \"logistic\".", type=str, default="logistic")
	parser.add_argument("--threads", help="Number of threads to use where applicable.", type=int, default=None)
	parser.add_argument("--subsets", help="Number of group-wise sub-sets to split the data into initially (Use 0 for automatic determination based on free memory).", type=int, default=None)
	parser.add_argument("--dropout", help="File containing list of feature indices to ignore.", type=str, default=None)
	parser.add_argument("--top_ft", help="Generate feature-level contribution heatmap using only the top N features.", type=int, default=None)
	parser.add_argument("--top_ft_window", help="Generate feature-level contribution heatmap using only the top 10 N feature windows.", type=int, default=None)
	parser.add_argument("--aim_max_iter", help="Maximum number of iterations to run in AIM mode.", type=int, default=10)
	parser.add_argument("--aim_max_ft", help="Maximum number of features to select in AIM mode.", type=int, default=1000)
	parser.add_argument("--aim_window", help="Window of top features to select from at each AIM iteration.", type=int, default=100)
	parser.add_argument("--aim_acc_cutoff", help="Minimum balanced accuracy a set of features must achieve to be selected.", type=float, default=0.9)
	parser.add_argument("--DrPhylo", help="Run Dr Phylo type analysis.", action='store_true', default=False)
	parser.add_argument("--AIM", help="Run ancestry informative markers analysis.", action='store_true', default=False)

	args = parser.parse_args()

	## Potentially load bearing cruft ##
	args.skip_preprocessing = False
	args.skip_processing = False
	args.preserve_xml = False
	args.fuzz_indels = False
	args.ensemble_parts = None
	args.ensemble_coverage = 5
	args.sparsify = False
	args.grid_threads = None
	## End cruft ##
	args.debug = True

	if args.threads is None:
		args.threads = os.cpu_count()
		if args.threads is not None:
			args.threads = math.floor(0.75 * args.threads)
		else:
			args.threads = 1
	args.threads = max(1, args.threads)
	if os.name == "posix":
		args.threads = 1  # Multithreading currently causes an issue when static compiling with glibc

	if args.classes is None and args.tree is None:
		raise Exception("Must invoke one of --tree or --classes option.")
	if args.classes is not None and args.tree is not None:
		raise Exception("Cannot use --tree and --classes options simultaneously.")
	if args.tree is None and args.gen_clade_list:
		raise Exception("Cannot use --gen_clade_list option and --classes options simultaneously.")
	#if not args.AIM and not args.DrPhylo:
	#	if args.lambda1_grid is None and args.lambda1 is None:
	#		args.lambda1 = 0.1
	if args.lambda1_grid is None:
		if args.lambda1 is not None:
			args.lambda1_grid = "{},{},{}".format(args.lambda1, args.lambda1 + 0.00001, 0.00002)
		elif args.AIM or args.DrPhylo:
			args.lambda1_grid = "0.1,1.0,0.1"
		else:
			# args.lambda1_grid = "0.1,0.10001,0.00002"
			# Default to 10x10 grid-search
			args.lambda1_grid = "0.1,1.0,0.1"
	if args.lambda2_grid is None:
		if args.lambda2 is not None:
			args.lambda2_grid = "{},{},{}".format(args.lambda2, args.lambda2 + 0.00001, 0.00002)
		else:
			if args.AIM:
				args.lambda2_grid = "0.0001,0.00011,0.0001"
			elif args.DrPhylo:
				args.lambda2_grid = "0.1,1.0,0.1"
			else:
				# args.lambda2_grid = "0.1,0.10001,0.00002"
				# Default to 10x10 grid-search
				args.lambda2_grid = "0.1,1.0,0.1"
	if args.DrPhylo:
		args.grid_summary_only = True
		if args.class_bal is None:
			if args.classes is not None:
				args.class_bal = "weighted"
			if args.tree is not None:
				args.class_bal = "phylo"
		if args.lambda1_grid is None:
			args.lambda1_grid = "0.1,1.0,0.1"
		if args.lambda2_grid is None:
			args.lambda2_grid = "0.1,1.0,0.1"
		if args.min_groups is None:
			args.min_groups = 1
		if args.stats_out == "":
			args.stats_out = "PGHS"
		if args.m_grid is None:
			args.m_grid = "20,30"
		if args.min_groups < 0:
			args.min_groups = 3
	if args.lambda1_grid is not None and args.lambda2_grid is not None and args.kfold > 1:
		raise Exception("Cannot use --kfold option while running in grid search mode.")
	args.single_lambda_pair = False
	if args.AIM:
		if args.lambda1_grid is None:
			args.lambda1_grid = "0.1,1.0,0.1"
		if args.lambda2_grid is None:
			args.lambda2_grid = "0.0001,0.00011,0.0001"
	#Convert single run to 1x1 grid run
	if args.lambda1_grid is None and args.lambda2_grid is None:
		args.lambda1_grid = "{},{},{}".format(args.lambda1, args.lambda1 + 0.00001, 0.00002)
		args.lambda2_grid = "{},{},{}".format(args.lambda2, args.lambda2 + 0.00001, 0.00002)
		# Disabled because it breaks subsetting
		# args.single_lambda_pair = True

	#sg_lasso still needs lambda1 and lambda2 values even when it ignores them for the lambda-pair file
	if args.lambda1 is None:
		args.lambda1 = float(args.lambda1_grid.split(",")[0])
	if args.lambda2 is None:
		args.lambda2 = float(args.lambda2_grid.split(",")[0])

	### Convert new parameter names to old ###
	args.partitions = args.subsets
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
		args.smart_sampling = 3
	elif args.class_bal == "phylo_1":
		args.smart_sampling = 1
	elif args.class_bal == "phylo_2":
		args.smart_sampling = 2
	args.grid_z = args.lambda1_grid
	args.grid_y = args.lambda2_grid
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
	args.timers = {"preprocessing": {"start": None}, "sglasso": {"start": None}, "analysis": {"start": None}}
	if (args.grid_z is None and args.grid_y is not None) or (args.grid_y is None and args.grid_z is not None):
		raise Exception("Only one grid search parameter specified, --lambda1_grid and --lambda2_grid must be specified together.")
	elif args.grid_z is not None and args.grid_y is not None and args.xval <= 1:
		gs_files = None
		if not os.path.exists(args.output):
			os.mkdir(args.output)
		args.aln_list = pf.aln_list_to_absolute(args.aln_list, args.output)
		if args.partitions is None or args.partitions == 1:
			if args.AIM and args.classes is not None:
				pf.aim_search(args)
			else:
				gs_files = pf.grid_search(args)
		elif args.partitions < 1:  # If partitions set to a value < 1, automatically determine ideal partition size
			try:  # run grid search on initial set, in a try block.
				gs_files = pf.grid_search(args)
			except Exception as e:  # if it fails, confirm that it failed because the dataset is too large for memory and set partitions accordingly, otherwise raise an error
				if str(e).split(":")[0] != "Minimum required subsets":
					raise Exception(e)
				else:
					args.partitions = math.ceil(float(str(e).split(":")[1]))
					args.disable_mc = True
					print("Data too large for available memory, setting subset count to {} and restarting.".format(args.partitions))
			# if it succeeds, you're done
		if args.partitions is not None and args.partitions > 1:
			if gs_files is None:
				gs_files = pf.generate_hypothesis_set(args)
			# fetch alignment filesizes to estimate split
			aln_stats = pf.stat_aln_files(args.aln_list)
			# generate partitions# alignment file list files and matching modified copies of args to match
			# aln_list_partitions = pf.generate_aln_partitions(aln_stats, args.partitions, args.aln_list)
			aln_list_partitions = pf.generate_aln_partitions(aln_stats, args.partitions, os.path.join(args.output, "aln_list.txt"))
			combo_count = int((args.partitions / 2) * (args.partitions - 1))
			# run grid_search separately for each hypothesis
			final_gs_files = [{} for j in range(0, len(gs_files['hypothesis_files']))]
			try:
				for j in range(0, len(gs_files['hypothesis_files'])):
					# run grid_search on each partition
					part_gs_files = [{} for i in range(0, combo_count)]
					part_args = copy.deepcopy(args)
					part_args.timers = args.timers
					part_args.stats_out = "G"
					part_args.grid_summary_only = True
					shutil.copy(gs_files['hypothesis_files'][j], os.path.join(args.output, os.path.basename(gs_files['hypothesis_files'][j]).replace("_hypothesis.txt".format(args.output), ".txt")))
					part_args.response = os.path.join(args.output, os.path.basename(gs_files['hypothesis_files'][j]).replace("_hypothesis.txt".format(args.output), ".txt"))
					part_args.kfold_ids = gs_files['hypothesis_files'][j]
					for i in range(0, combo_count):
						part_args.output = os.path.join(args.output, "{}_part{}".format(os.path.basename(args.output), i))
						part_args.partitions = 1
						part_args.aln_list = aln_list_partitions[i]
						part_gs_files[i] = pf.grid_search(part_args)
						# if i == args.partitions + 1:
						# 	os.remove(aln_list_partitions[i])
					# analyze GSS file out of grid_search to select genes that capture X% of total GSS
					selected_size = sum(aln_stats.values())
					selected_proteins = []
					selected_groups = []
					pct_select = 100.0
					while selected_size > (sum(aln_stats.values()) / args.partitions):
						pct_select = pct_select * 0.9
						selected_proteins = []
						selected_groups = []
						for i in range(0, combo_count):
							GSS_file = part_gs_files[i]['GSS_median_files'][0]
							selected_proteins += pf.select_top_gss(GSS_file, pct_select)
						selected_size = 0
						for file_group in aln_stats.keys():
							for filename in file_group.strip().split(','):
								if os.path.splitext(os.path.basename(filename))[0] in selected_proteins:
									selected_groups += [file_group]
									selected_size += aln_stats[file_group]
									break
					# generate final alignment file list file and matching modified copy of args to match
					aln_dir = os.path.split(args.aln_list)[0]
					final_aln_list_filename = os.path.join(aln_dir, "{}_aln_list_final{}".format(os.path.splitext(os.path.basename(part_args.response))[0], ".txt"))
					with open(final_aln_list_filename, 'w') as file:
						for group in selected_groups:
							file.write("{}\n".format(group))
					# run grid search on final set.
					part_args.output = os.path.join(args.output, "{}_final".format(args.output))
					part_args.aln_list = final_aln_list_filename
					part_args.disable_mc = True
					part_args.grid_summary_only = args.grid_summary_only
					part_args.stats_out = args.stats_out
					final_gs_files[j] = pf.grid_search(part_args)
					os.remove(gs_files['hypothesis_files'][j].replace("_hypothesis.txt".format(args.output), ".txt"))
					for i in range(0, combo_count):
						pf.cleanup_directory(args, part_gs_files[i])
					pf.cleanup_directory(args, final_gs_files[j])
				# pf.cleanup_directory(args, gs_files)
				for filename in os.listdir(args.output):
					if os.path.splitext(filename)[1] == ".txt" and os.path.exists(os.path.join(args.output, filename)) and not args.debug:
						os.remove(os.path.join(args.output, filename))
			except Exception as e:
				raise Exception(e)
			# finally:
			# 	# clean up partition folders
			# 	for i in range(0, args.partitions):
			# 		try:
			# 			shutil.move("{}_part{}".format(args.output, i), os.path.join(args.output, "{}_part{}".format(args.output, i)))
			# 		except:
			# 			pass
			# 	try:
			# 		shutil.move("{}_final".format(args.output), os.path.join(args.output, "{}_final".format(args.output)))
			# 	except:
			# 		pass
	for key in args.timers.keys():
		print("Total time elapsed for {}: {}.".format(key, args.timers[key]["total"]))
