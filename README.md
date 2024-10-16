# Evolutionary Sparse Learning with MyESL #

## Introduction ##
MyESL is a set of programs written in C++ and Python designed to build an evolutionary sparse learning model for a phylogenetic hypothesis or binary classes of species with or without a phenotype. MyESL uses sparse group lasso regression to optimize weights for each position in the sequence alignments, quantifying the strength of association between a position and the hypothesis/class used in the analysis. MyESL also pre-posses input data and processes ESL model outputs for further biological discoveries. MyESL is currently available as a standalone software `(MyESL.exe)` compatible with Windows and Linux operating systems. Moreover, it offers the option of installation via a Python environment to conduct ESL analysis. 

## Implementation of MyESL ## 
We discuss every required and optional argument for ESL analysis of an example dataset using MyESL software (MyESL.exe). Next, we will discuss installing MyESL in the Python environment.  

### Download ###
You can download `MyESL` from the GitHub repository or use the command  

	git clone https://github.com/kumarlabgit/MyESL MyESL
	cd MyESL-master

 
### Required arguments:

<br />

```

alignment_list.txt                       : A text file contains a list of paths for all sequence alignments. For example, sequence alignments from the fungi dataset provided within the `Fungi_data` directory:
                                           aln\BUSCOfEOG7008JB.fasta 
                                           aln\BUSCOfEOG7008JM.fasta 
                                           aln\BUSCOfEOG7008JT.fasta 

--tree <phylogenetic_tree.nwk>           : A phylogenetic tree in newick format with a node ID to construct a hypothesis for the clade of interest.
                                           The hypothesis can also be specified with a separate file using the --classes parameter.
                                           It is highly recommended that the number of species in the clade be equal to or less than those outside of the clade.
                                           It is also recommended to use the smart sampling option (—-class_bal) when the number of species inside the clade is less or greater than the number outside the clade.
OR

--classes <classes.txt>                  : Requires a text file containing a user-defined hypothesis. It has two columns, which are tab-separated. The first column contains species names, and the second column contains the response 
                                          value for the species (+1/-1). A member species in the clade or with a specific phenotype receives +1 and -1 otherwise. This hypothesis is unconstrained by the tree structure. It is highly 
                                          recommended that the number of species with the response +1 equal the number of species with -1. 

```
<br />	

### Optional arguments:

Users can also specify other options in MyESL for processing the input data, building ESL models, and post-processing the ESL models. 

#### Pre-processing input data:

<br />

```
--clade_list <filename.txt>         : Users can test multiple phylogenetic hypotheses when the input phylogenetic tree contains multiple clade IDs. Clade list file should have one clade ID on each line. This option must be used with "--tree" option.

--gen_clade_list <int,int>          : Users can generate multiple hypotheses when the input phylogeny contains no clade ID. The size of the clade is determined by the input integers defining the upper and lower limits
                                      of clade size, respectively.
   
--class_bal <string>                : MyESL performs class balancing, a common practice in supervised machine learning. Class balancing helps balance the number of species inside and outside the focal clade of interest.
                                      Class balancing in MyESL is performed by phylogenetic aware sampling <phylo> when the "--tree" option provides the phylogenetic hypothesis. MyESL also makes a balance between classes 
                                      by using inverse weights using the option <weight>. Upsampling and downsampling between classes can also be performed in MyESL using <up> and <down> options. 

--data_type <string>                : <nucleotide> informs MyESL to treat A, T, C, G, and U as valid characters without case sensitivity, and all other characters will be treated as missing data.
                                      <protein> option treats all unambiguous IUPAC amino acid letters (case insensitive) as valid characters.
                                      <molecular> option provides a way to use both nucleotide and acid letters as valid characters.
                                      <universal> option is used as default, which allows the analysis of presence/absence (0/1) data.

--bit_ct <int>                      : One can drop all bit-columns in which the bit 1 appears fewer than a certain number of times.

```
<br />	

#### ESL model building:

<br />

```
--lambda1 <float>                    : The site sparsity parameter that ranges from 0 to 1. When not specified, the default is 0.1. It is required to build a single-clade model.

--lambda2 <float>                    : The gene sparsity parameter ranges from 0 to 1, and the default is 0.1 when not specified. It is required to build a single-clade model.

--lambda1_range <min,max,step>        : This option allows users to set the range for the site sparsity parameter. The site sparsity grid is defined by a string of float numbers min, max, and step_size, which range from 0 to 1.
                                       For example, --lambda1_grid 0.1,0.9,0.1. This option must be used with --lambda2_grid.  

--lambda2_grid <min,max,step>        : This option allows users to set the range for the group sparsity parameter. The group sparsity grid is defined by a string of float numbers min, max, and step_size, which range from 0 to 1.
                                       For example, --lambda2_grid 0.1,0.9,0.1. This option must be used with --lambda1_grid. 

--min_groups <int>                   : This option allows users to set the minimum number of genes included in the multi-gene ESL models and helps early stopping in the grid search over the sparsity parameter space.
                                       It takes a value greater than zero (0) and builds models containing more or equal numbers of groups in the model.

--group_wt <filename.txt>            : MyESL uses the square root of the number of bit columns for groups as group weights by default. Users can provide group weights via a text file (<filename.txt>), a tab-separated two-column file.
                                       The first column contains group names, and the second contains corresponding group weights.

--kfold <int>                        : This directive allows users to build multiple ESL models using k-fold cross-validation. For example, if k-fold is set to 5, 80% of rows are used in model training, and 20% of rows are withheld for validation.

--no_group_penalty                   : This directive allows users to perform LASSO regression without group penalty. The LASSO regression without group penalty is also performed when a single FastA file is provided as input. 

```
<br />	

#### Processing ESL-model outputs:

<br />

```
--output <string>                   : The name of the output directory where all results from MyESL analysis will be stored. The program creates an output directory automatically if it is not specified.

--stats_out <PGHS>                  : MyESL processes the ESL model and outputs, Position Sparsity Scores (<P>), Gene Sparsity Scores (<G>), Hypothesis Sparsity Scores (<H>), and
                                      Species Sparsity Scores and Prediction Probability (<S>). Users can output multiple files using multiple inputs like <BPS>. Details of these scores can be found in reference #1.  
```
<br />	

### MyESL output files

MyESL produces multiple output files based on different user analysis directives. 

#### Output files from a single ESL model

<br />

```
ESL_model_{clade_ID}.txt   : A tab-separated text file containing allele at each position in a group and non-zero beta values estimated from sparse group lasso analysis. 

PSS_{clade_ID}.txt         : A tab-separated text file containing position ID in a group and non-zero position sparsity scores.
GSS_{clade_ID}.txt         : A tab-separated text file containing group ID and non-zero group sparsity scores.
HSS_{clade_ID}.txt         : A text file containing hypothesis sparsity scores.
SPS_SPP_{clade_ID}.txt     : A tab-separated text file containing species name, regression response, species prediction score, and prediction probability. 

```
<br />

#### Output files from multiple ESL models

When grid options for penalty parameters are used, MyESL produces multiple ESL models and other sparsity scores. In addition to a single model output for each lambda pair, MyESL outputs a summary of all single models. 

<br />

```
ESL_model_{clade_ID}_l1_l2.txt   : A tab-separated text file containing allele at each position in a group and non-zero beta values estimated from sparse group lasso analysis. 

PSS_{clade_ID}_l1_l2.txt         : A tab-separated text file containing position ID in a group and non-zero position sparsity scores.
GSS_{clade_ID}_l1_l2.txt         : A tab-separated text file containing group ID and non-zero group sparsity scores.
HSS_{clade_ID}_l1_l2.txt         : A text file containing hypothesis sparsity scores.
SPS_SPP_{clade_ID}_l1_l2.txt     : A tab-separated text file containing species name, regression response, species prediction score, and prediction probability. 

PSS_{clade_ID}_summary.txt         : A tab-separated text file containing position ID in a group and non-zero position sparsity scores.
GSS_{clade_ID}_summary.txt         : A tab-separated text file containing group ID and non-zero group sparsity scores.
HSS_{clade_ID}_summary.txt         : A text file containing hypothesis sparsity scores.
SPS_SPP_{clade_ID}_summary.txt     : A tab-separated text file containing species name, regression response, species prediction score, and prediction probability. 

Note: l1 and l2 refer to site and group sparsity parameters used for model building in grid search.
```
<br />

### Usage ###

After downloading and setting the current working directory to MyESL, one can perform ESL analysis using `MyESL.exe.` The Fungi_data directory provides an example dataset. The directory contains 1,232 gene sequence alignments in fasta format. It also contains a phylogenetic tree of 86 fungi with clade IDs for two clades. The clade ID `Clade_X1` indicates a clade found to be fragile in a previous study (Sharma and Kumar 2024), and the clade `Control` was analyzed as a control clade previously.

<br />

```
MyESL.exe Fungi_data\aln.txt  --tree Fungi_data\Fungi_T1_with_ID.nwk --lambda1 0.1 --lambda2 0.2 --output Fungi_out

OR

MyESL.exe alignment_list.txt  --classes classes.txt
```
<br />

This command will produce two models for clades "Clade_X1" and "Control". The ESL model for individual clades can also be built using the option `--clade_list`. The clade model for the clade `Clade_X1` and `Control` can be made using the following command, respectively:

```
MyESL.exe Fungi_data\aln.txt  --tree Fungi_data\Fungi_T1_with_ID.nwk --clade_list\clade_X1.txt --lambda1 0.1 --lambda2 0.2 --output Fungi_out_clade_X1
MyESL.exe Fungi_data\aln.txt  --tree Fungi_data\Fungi_T1_with_ID.nwk --clade_list\clade_Control.txt --lambda1 0.1 --lambda2 0.2 --output Fungi_out_clade_Control

```

It is always recommended to keep the number of species in both classes inside and outside the clade the same or at least approximately similar. However, if the number of species within the clade differs from those outside, MyESL provides a way to balance the species count within the clade. Two approaches are available for this purpose: phylogenetically aware subsampling for class balancing and class weighting. When the number of species within the clade exceeds that outside, class weighting is the preferred method. 


Phylogenetic-aware class balancing:

```
MyESL.exe Fungi_data\aln.txt  --tree Fungi_data\Fungi_T1_with_ID.nwk --clade_list\clade_X1.txt --lambda1 0.1 --lambda2 0.2 --class_bal phylo --output Fungi_out_clade_X1
```
Class balancing using inverse class weights:

```
MyESL.exe Fungi_data\aln.txt  --tree Fungi_data\Fungi_T1_with_ID.nwk --clade_list\clade_X1.txt --lambda1 0.1 --lambda2 0.2 --class_bal weight --output Fungi_out_clade_X1
```

An ESL model can also be built using a text-based response file without requiring a phylogenetic tree in Newick format. In the response file, species within a clade or those with a predefined trait are assigned a value of 1, while all others are assigned -1. A response file was created for the clade "Clade_X1" in the `Fungi_data` directory. 

```
MyESL.exe Fungi_data\aln.txt  --classes Fungi_data\Clade_X1.txt --clade_list\clade_Control.txt --lambda1 0.1 --lambda2 0.2 --class_bal weight --output Fungi_out_clade_X1
```

In the previous examples, all models were built using predefined site and group sparsity parameter values. Users can also create a consensus ESL model by combining multiple models built with sparsity parameter values ranging from 0 to 1, with a step size of 0.1. Additionally, users can customize the range of lambda values and the step size.

An ensemble ESL model without specifying the range of sparsity parameter values:

```
MyESL.exe Fungi_data\aln.txt  --tree Fungi_data\Fungi_T1_with_ID.nwk --clade_list\clade_X1.txt --class_bal weight --output Fungi_out_clade_X1
```

An ensemble ESL model by specifying the range of sparsity parameter values. The site (`lambda1`) and group (`lambda2`) sparsity parameters with a range from 0.1 to 0.9 with a step size of 0.1:

```
MyESL.exe Fungi_data\aln.txt  --tree Fungi_data\Fungi_T1_with_ID.nwk --lambda1_grid 0.1,0.9,0.1 --lambda2_grid 0.1,0.9,0.1 --clade_list\clade_X1.txt --class_bal weight --output Fungi_out_clade_X1
```

MyESL produces multiple numeric outputs and a graphical output. The numeric outputs include different sparsity scores that quantify the association between features (e.g., bits, sites, genes) and hypotheses tested. The model also calculates the hypothesis sparsity scores (HSS), providing overall support for the hypothesis tested. MyESL usually produces all these scores in a text file format by default. However, users can specify which sparsity scores need to be produced by MyESL using the option `--stats.`

MyESL outputs only gene and site sparsity scores by using the argument `--stats GS`:

```
MyESL.exe Fungi_data\aln.txt  --tree Fungi_data\Fungi_T1_with_ID.nwk --clade_list\clade_X1.txt --lambda1 0.1 --lambda2 0.2 --class_bal phylo --stats SG --output Fungi_out_clade_X1
```

The graphical output is a grid, with each row representing the species (with normalized classification probability within the parenthesis) used in the analysis and columns for genes selected by the ESL model. Each cell in the grid is the gene-species concordance (GSC) score. A positive value for the gsc implies that the gene supports the placement of the species within the clade and is presented in green. In contrast, the negative value indicates discordance, represented by red. Users can modify the grid output by specifying the number of genes and the number of species to be displayed using the option ``--m_grid``. In this case, only a specified number of species from the clade of interest will be displayed and sorted using the classification probability. 

```
MyESL.exe Fungi_data\aln.txt  --tree Fungi_data\Fungi_T1_with_ID.nwk --clade_list\clade_X1.txt --lambda1 0.1 --lambda2 0.2 --class_bal phylo --stats SG --m_grid 20,20 --output Fungi_out_clade_X1
```

### DrPhylo analysis using MyESL

Users can perform DrPhylo analysis using MyESL. Instructions and Outputs for DrPhylo analysis will be found at [DrPhylo](https://github.com/kumarlabgit/MyESL/tree/DrPhylo). 

### Installation of MyESL into Python ##

To run MyESL, you will need Python 3.8 or later installed, as well as the following Python libraries:
```R
numpy
biopython
matplotlib
pandas
```

You can install these libraries using pip:

`pip install biopython numpy pandas matplotlib`

To perform DrPhylo analysis, get the ESL package using the following on the command line:

	git clone -b master https://github.com/kumarlabgit/MyESL MyESL
	cd MyESL
	bash setup.sh

#### MyESL pipeline uses the same directives as MyESL.exe. To perform an analysis using the example dataset, only replacing `MyESL.exe` with `MyESL.py`

## References ##
If you use MyESL in your research, please cite our articles:

1. Kumar, S. and Sharma, S (2021). Evolutionary Sparse Learning for Phylogenomics, Molecular Biology and Evolution, Volume 38, Issue 11, November 2021, Pages 4674–4682.
2. Sharma, S. & Kumar, S. (2024). Discovering fragile clades and causal sequences in phylogenomics by evolutionary sparse learning. (In review)
3. Sanderford et al., (2024).  MyESL: A software for evolutionary sparse learning in molecular phylogenetics and genomics (In preparation).
