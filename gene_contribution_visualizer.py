import os
import argparse
import shutil
import math
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import TwoSlopeNorm
import numpy as np

def get_roc(labels, predictions):
	span = float(max(labels) - min(labels))
	pos_set = [x[1] for x in zip(labels, predictions) if x[0] == max(labels)]
	neg_set = [x[1] for x in zip(labels, predictions) if x[0] == min(labels)]
	TPR = []
	FPR = []
	for cutoff in [(x * (span / 100)) + min(labels) for x in range(1, 99)]:
		TPR.append(sum([1.0 for x in pos_set if x >= cutoff]) / len(pos_set))
		FPR.append(sum([1.0 for x in neg_set if x >= cutoff]) / len(neg_set))
	print("Accuracy at cutoff = 0:")
	print("True Positive Rate: {}".format(sum([1.0 for x in pos_set if x >= 0]) / len(pos_set)))
	print("True Negative Rate: {}".format(sum([1.0 for x in neg_set if x < 0]) / len(neg_set)))
	print("False Positive Rate: {}".format(sum([1.0 for x in neg_set if x >= 0]) / len(neg_set)))
	print("False Negative Rate: {}".format(sum([1.0 for x in pos_set if x < 0]) / len(pos_set)))
	return TPR, FPR


def merge_predictions(file_list, output_filename=None):
	seqid_list = None
	data = {}
	headers = {}
	fields = set()
	for predictions_table in file_list:
		with open(predictions_table, 'r') as file:
			headers[predictions_table] = [val.strip() for val in file.readline().split("\t")]
			for field in headers[predictions_table]:
				fields.add(field)
			file.seek(0)
			data[predictions_table] = np.genfromtxt(file, dtype=None, skip_header=1, missing_values="N/A", encoding=None)
			if seqid_list is not None:
				if not (seqid_list==np.array([[val[0] for val in data[predictions_table]]]).T).all():
					raise Exception("Cannot combine prediction tables for differing sequence sets.")
			else:
				seqid_list = np.array([[val[0] for val in data[predictions_table]]]).T
	output_header = copy.copy(list(headers.values())[0])
	for field in list(fields):
		if field not in output_header:
			output_header.append(field)
	for predictions_table in file_list:
		data[predictions_table] = np.asarray([list(val)[1:] for val in data[predictions_table]], dtype="float")
#		print(headers[predictions_table])
		headers[predictions_table] = headers[predictions_table][1:]
#		print(headers[predictions_table])
	for field in output_header:
		if field=="SeqID":
			#merged_data = seqid_list
			col_list = [seqid_list]
		else:
#			print(field)
#			print([headers[predictions_table].index(field) for predictions_table in file_list if field in headers[predictions_table]])
			cols = np.array([data[predictions_table][:,headers[predictions_table].index(field)] for predictions_table in file_list if field in headers[predictions_table]])
			col = np.array([np.mean(cols, 0)]).T
			col_list.append(col)
#			print(col)
#			print(merged_data)
#			merged_data = np.concatenate(merged_data, col)
#	print(col_list)
#	for pair in zip(output_header, col_list):
#		print("Field:{}\tLength:{}\tShape:{}".format(pair[0],len(pair[1]), pair[1].shape))
#	merged_data = np.rec.fromarrays(col_list, names=",".join(output_header))
	merged_data = np.hstack(tuple(col_list))
	if output_filename is None:
		output_filename="testing.txt"
#	print(len(output_header))
#	print(merged_data.dtype)
#	print(merged_data.shape)
#	np.savetxt(output_filename, merged_data, delimiter="\t", fmt='%s,'+','.join(['%.2f']*(merged_data.shape[1]-1)))
	np.savetxt(output_filename, merged_data, fmt='\t'.join(['%s']*merged_data.shape[1]), header='\t'.join(output_header), comments='')
#	np.savetxt(output_filename, merged_data, delimiter="\t")
	return output_filename


def main(predictions_table, lead_cols=4, response_idx=2, prediction_idx=3, output=None, ssq_threshold=0, gene_limit=100, m_grid=False):
	sample_size = 24
	sample_size = None
	with open(predictions_table, 'r') as file:
		header = file.readline().split("\t")
		num_cols = len(header)
		file.seek(0)
		if sample_size is None:
			data = np.genfromtxt(file, dtype='U', skip_header=1, usemask=True, missing_values="N/A", encoding=None)
		else:
			data = np.genfromtxt(file, dtype='U', skip_header=1, usemask=True, missing_values="N/A", encoding=None, usecols=range(0, lead_cols + sample_size))
			header = header[0: lead_cols + sample_size]
		seqid_list = [(val)[0] for val in data]
		data = np.asarray([list(val)[1:] for val in data], dtype="float")
	num_rows, num_cols = data.shape
	print("\nRows: {}\nCols: {}\n".format(num_rows, num_cols))


	# Sort columns by sum of contributions
	ssq_scores = list(zip(range(lead_cols-1, num_cols), [np.sum(np.nan_to_num(gene_col)**2) for gene_col in data.transpose()[lead_cols-1:]], header[lead_cols:]))
	ssq_scores.sort(key=lambda tup: tup[1], reverse=True)
	if ssq_threshold > 0:
		ssq_scores = [val for val in ssq_scores if val[1] >= ssq_threshold]
	sum_scores = list(zip(range(lead_cols-1, num_cols), [np.sum(gene_col) for gene_col in data.transpose()[lead_cols-1:]], header[lead_cols:]))
	sum_scores.sort(key=lambda tup: tup[1], reverse=True)
	ingroup_data = data[data[:, 0] == 1.0]
	# Mean positive GSC scores
	mpgsc_scores = list(zip(range(lead_cols-1, num_cols), [np.mean([x for x in np.nan_to_num(gene_col) if x > 0] if len([x for x in np.nan_to_num(gene_col) if x > 0]) > 0 else [0]) for gene_col in ingroup_data.transpose()[lead_cols-1:]], header[lead_cols:]))
	mpgsc_scores.sort(key=lambda tup: tup[1], reverse=True)


	if m_grid:
		data = data[:, list(range(0, lead_cols - 1)) + [val[0] for val in mpgsc_scores[0:gene_limit]]]
		header = [val[2] for val in mpgsc_scores[0:gene_limit]]
	else:
		data = data[:, list(range(0, lead_cols - 1)) + [val[0] for val in ssq_scores[0:gene_limit]]]
		header = header[1:lead_cols] + [val[2] for val in ssq_scores[0:gene_limit]]

	# Sort rows by expit of predicted value
	expit = [1/(1 + math.exp(-x)) for x in data[:, 1]]
	pr = [(x - 0.5)/(max(expit) - 0.5) for x in expit]
	pr = [x if x > 0 else 0 for x in pr]
	sorted_rows = list(zip(range(0, num_rows), pr, data[:, 0]))
	sorted_rows.sort(key=lambda tup: (tup[2], tup[1]), reverse=True)
	data = data[[val[0] for val in sorted_rows]]
	seqid_list = ["{} ({:0.2f})".format((val)[0], val[1]) for val in zip(seqid_list, pr)]

	# temp1 = list(range(lead_cols-1, len(data[0])))
	# print(temp1)
	# print(len(temp1))
	# print(len(data))
	# print(len(data[0]))

	if m_grid:
		data = data[0:sum([1 for val in data[:, 0] if float(val)>0.99]), list(range(lead_cols-1, len(data[0])))]

	seqid_list = [seqid_list[val[0]] for val in sorted_rows[0:len(data)]]

	# Reset dimensions if truncation has occurred.
	num_rows, num_cols = data.shape

	# draw gridlines
	xtick_width = 20
	ytick_width = 20
	fig, ax = plt.subplots(figsize=(num_cols * 0.1, num_rows * 0.1))
	label_size = 3
	cell_label_size = 2.0
	DPI = 600
	ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=0.2)
	# Make red-yellow-green color map with NaN=blue
	cmap = copy.copy(mpl.cm.get_cmap("RdYlGn"))
	cmap.set_bad(color='lightskyblue')
	# Apply color map to first three columns, and then separately to the rest of the columns
	if m_grid:
		ax.imshow(data[:, :], cmap=cmap, norm=TwoSlopeNorm(0), extent=(0, num_cols * xtick_width, num_rows * ytick_width, 0))
	else:
		# ax.imshow(data[:, 0:3], cmap=cmap, norm=TwoSlopeNorm(0), extent=(0, (lead_cols - 1) * xtick_width, num_rows * ytick_width, 0))
		# ax.imshow(data[:, 3:], cmap=cmap, norm=TwoSlopeNorm(0), extent=((lead_cols - 1)*xtick_width, num_cols*xtick_width, num_rows*ytick_width, 0))
		#for i in range(0, lead_cols-1):
		for i in range(0, 2):
			ax.imshow(data[:, i:i+1], cmap=cmap, norm=TwoSlopeNorm(0), extent=(i*xtick_width, (i+1) * xtick_width, num_rows * ytick_width, 0))
		ax.imshow(data[:, lead_cols-1:], cmap=cmap, norm=TwoSlopeNorm(0), extent=((lead_cols - 1)*xtick_width, num_cols * xtick_width, num_rows * ytick_width, 0))
	ax.set_xticks(np.arange(0, num_cols * xtick_width, xtick_width));
	ax.set_xticklabels(header, rotation=90, ha='left', size=label_size)
	ax.set_yticks(np.arange(0, num_rows * ytick_width, ytick_width));
	ax.set_yticklabels(seqid_list, va="top", size=label_size)
	ax.set_xlabel('Alignment Names', fontsize=label_size*1.25)
	ax.set_ylabel('Sequence IDs', fontsize=label_size*1.25)
	ax.set_title('Group Contribution Weights', fontsize=label_size*1.75)
	for (i, j), z in np.ndenumerate(data):
		if not m_grid:
			if j < lead_cols - 1:
				ax.text((j+0.5) * xtick_width, (i+0.5)*ytick_width, '{:0.1f}'.format(z), ha='center', va='center', size=cell_label_size)
			else:
				ax.text((j + 0.5) * xtick_width, (i + 0.5) * ytick_width, '{:0.2f}'.format(z), ha='center', va='center', size=cell_label_size)
		elif np.isnan(z):
			ax.text((j + 0.5) * xtick_width, (i + 0.5) * ytick_width, 'x', ha='center', va='center', size=cell_label_size)
	if output is None:
		output = os.path.join(os.path.dirname(predictions_table),"{}.png".format(os.path.splitext(os.path.basename(predictions_table))[0]))
	plt.savefig(output, dpi=DPI, bbox_inches='tight')
	plt.close()
	try:
		roc = get_roc(data[:,0], data[:,1])
		roc_output = os.path.join(os.path.dirname(predictions_table),"{}_ROC.png".format(os.path.splitext(os.path.basename(predictions_table))[0]))
		plt.plot(roc[1], roc[0])
		plt.savefig(roc_output, dpi=DPI, bbox_inches='tight')
	except:
		print("Couldn't generate ROC output for {}.".format(predictions_table))
	return output


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Gene contribution visualizer for ESL pipeline.")
	parser.add_argument("predictions_table", help="Table of gene prediction contribution values to sort and color.", type=str)
	parser.add_argument("--lead_cols", help="Number of columns ahead of gene contribution columns.", type=int, default=4)
	parser.add_argument("--response_idx", help="1-based index of response column.", type=int, default=2)
	parser.add_argument("--prediction_idx", help="1-based index of prediction column.", type=int, default=3)
	parser.add_argument("--output", help="Output image file.", type=str, default=None)
	parser.add_argument("--threshold", help="Drop columns with sum of squares below threshold value.", type=float, default=0)
	parser.add_argument("--m_grid", help="Make DrPhylo output.", action='store_true', default=False)
	args = parser.parse_args()
	main(args.predictions_table, args.lead_cols, args.response_idx, args.prediction_idx, args.output, args.threshold, m_grid=args.m_grid)
