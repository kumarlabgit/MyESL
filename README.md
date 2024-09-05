# DrPhylo analysis for detecting fragile clades in an inferred phylogeny from phylogenomic datasets #

DrPhylo investigates major species relationships in phylogenies inferred from phylogenomic datasets. This approach applies evolutionary sparse learning (ESL) to build a genetic model for a clade of interest for a given set of sequence alignments from genes or genomic loci (see citation #2). DrPhylo outputs:
<br />


```
Gene species concordance (GSC)              : A numeric value that identifies gene-species combinations (sequence) harboring concordant ( GSC > 0) or conflicting (GSC < 0) phylogenetic signals. 

Sequence Classification Probability (SCP)   : Classification probability for a species for being a member of the clade of interest.  

Clade Proabaility (CP)                      : Probability (0-1) for the clade of interest. A low value for CP indicates the plausibility of being fragile.  

Model grid                                  : This is a grid representation of the clade model. The rows of this grid represent species with SCP and columns for genes. 
                                              Green       (GSC > 0 )
                                              Red         (GSC < 0 )
                                              White       (GSC = 0 )
                                              Cross mark  (Missing data)    

```
<br />
A schematic outline of DrPhylo analysis. The figure is used from DrPhylo article (see citation #1) :
			<div style="display: flex; justify-content: center;">
			    <img src="https://github.com/ssharma2712/DrPhylo/assets/11808951/332cbd52-a1b3-4593-a0dd-62f6723376a0" width="600">
			</div>

## DrPhylo Analysis ## 

DrPhylo is integrated into MyESL, which can be implemented using the ```--DrPhylo``` flag in MyESL. Therefore, installing MyESL will allow users to implement DrPhylo. MyESL can be installed in the Python environment using the source codes and the instructions below. A user can also use MyESL as a standalone executable software, which does not require installation.    

## MyESL.exe for DrPhylo analysis ##

To run DrPhylo integrated into MyESL, download MyESL using: 

```
git clone -b DrPhylo https://github.com/kumarlabgit/MyESL DrPhylo
cd MyESL

```

## Basic input for DrPhylo analysis using default options ##
DrPhylo analysis in MyESL requires the list of paths for all groups (sequence alignments of genes or genetic loci in fasta format) as a text file and a phylogenetic hypothesis. The phylogenetic hypothesis tested by DrPhylo can be provided using a phylogenetic tree in newick format or a text file containing the hypothesis. 

## Implementation ##

After downloading, a user can perform DrPhylo analysis using the required inputs and other necessary optional arguments.

<br />

```
MyESL.exe alignment_list.txt  --tree phylogenetic_tree.nwk --DrPhylo

OR

MyESL.exe alignment_list.txt  --classes phylogenetic_hypothesis.txt --DrPhylo
```
<br />


#### Required arguments:

<br />

```

alignment_list.txt                       : A text file contains a list of paths for all sequence alignments. For example,
                                           angiosperm_alns/7276_C12.fasta
                                           angiosperm_alns/5111_C12.fasta
                                           angiosperm_alns/5507_C12.fasta

--tree phylogenetic_tree.nwk             : A phylogenetic tree in newick format with a node ID to construct a hypothesis for the clade of interest.
                                           The hypothesis can also be specified with a separate file using the --classes parameter.
                                           It is highly recommended that the number of species in the clade be equal to or greater than those outside of the clade.
                                           It is also recommended to use the smart sampling option (—-balancing) when the number of species inside the clade is greater than the number outside the clade.
OR

--classes phylogenetic_hypothesis.txt   : Requires a text file containing a user-defined hypothesis. It has two columns, which are tab-separated. The first column contains species names, and the second column contains the response value for the 
                                            species (+1/-1). A member species in the clade receives +1 and -1 otherwise. This hypothesis is unconstrained by the tree structure. It is highly recommended that the number of species within the clade of 
                                            interest (+1) is equal to the number of species outside the clade. The hypothesis can also be specified using a separate text file provided using the --response parameter.  

```
<br />	


DrPhylo builds multiple ESL models by performing a ```grid search``` over the discrete sparsity parameters (group and site) space. DrPhylo achieves computational efficiency by early terminating the grid-search process and selecting only multi-gene models. 

DrPhylo outputs a model grid (```M-grid```) and a text file in a matrix format containing the model. A ```M-grid``` is a two-dimensional graphical presentation of the ESL model containing species names with classification probability (rows) and groups sorted by influence (columns).  

#### Optional argumnets:

<br />

```

--lamda1_grid                : This option allows users to set the range for the site sparsity parameter. The site sparsity grid is defined by a string of float numbers min, max, step_size which range from 0 to 1.
                               For example, --lamda1_range 0.1, 0.9, 0.1. This option must be used with --lamda2_range.  

--lamda2_grid                : This option allows users to set the range for the group sparsity parameter. The group sparsity grid is defined by a string of float numbers min, max, step_size which range from 0 to 1.
                               For example, --lamda2_range 0.1, 0.9, 0.1. This option must be used with --lamda1_range. 

--min_groups                 : This option allows users to set the minimum number of genes included in the multi-gene ESL models and helps early stopping in the grid search over the sparsity parameter space.
                               It takes a value greater than zero (0) and builds models containing more or equal numbers of groups in the model.

--class_bal                  : DrPhylo also performs class balancing, a common practice in classification analysis of supervised machine learning. Class balancing helps to make a balance between the number of species inside and outside the focal clade of 
                               interest. Class balancing in DrPhylo is performed by phylogenetic aware sampling <phylo> when the phylogenetic hypothesis is provided by the "--tree" option. DrPhylo also makes a balance between classes by using inverse weights using
                               the option <weight>. 
--output                     : The name of the output directory where all results from DrPhylo analysis will be stored. The program creates this directory automatically. 
```
<br />	

## DrPhylo implementation using an example Fungi dataset in Windows ##

After finishing the setup, change the directory to `DrPhylo-master`. To create a text file containing the list of paths for all gene sequence alignments: 
```
cd Fungi_data
for %f in (aln\*.fasta) do echo %f >> aln.txt
cd ..
```
<br />
<br />

DrPhylo analysis for building clade models for two labeled clades, producing a result set for each of two hypotheses generated from the input tree:
<br />
<br />

```
MyESL.exe Fungi_data/aln.txt --tree Fungi_data/Fungi_T1.nwk --DrPhylo --output Fungi_tree_output  

```
<br />

DrPhylo analysis using a user-defined hypothesis, producing a single clade model for the given hypothesis:
<br />

```
MyESL.exe Fungi_data/aln.txt --classes Fungi_data/A_B_Hyp.txt --DrPhylo --output Fungi_test_output 

```
<br />

DrPhylo analysis for building an ensemble clade model for a user-defined hypothesis using the grid search option. The site and gene sparsity score range from 0.05 to 0.1 and a step size of 0.05:
<br />

```
MyESL.exe Fungi_data/aln.txt --classes Fungi_data/A_B_Hyp.txt --DrPhylo --lambda1_grid 0.05,0.1,0.05 --lambda2_grid 0.05,0.1,0.05 --output Fungi_test_grid_output 

```
<br />

## DrPhylo implementation using an example Fungi dataset in Windows ##

After finishing the setup, change the directory to `DrPhylo-master`. To create a text file containing the list of paths for all gene sequence alignments: 
```
cd Fungi_data
ls angiosperm_alns/*.fasta > angiosperm_100_sample_alns.txt
cd ..
```
<br />
<br />

DrPhylo analysis for building clade models for two labeled clades, producing a result set for each of two hypotheses generated from the input tree:
<br />
<br />

```
MyESL.exe sample_files/angiosperm_100_sample_alns.txt --tree sample_files/ESL_test.nwk --DrPhylo --output sample_tree_output  

```
<br />

DrPhylo analysis using a user-defined hypothesis, producing a single clade model for the given hypothesis:
<br />

```
MyESL.exe sample_files/angiosperm_100_sample_alns.txt --classes sample_files/test_pred.txt --DrPhylo --output sample_tree_output 

```
<br />

DrPhylo analysis for building an ensemble clade model for a user-defined hypothesis using the grid search option. The site and gene sparsity score range from 0.05 to 0.1 and a step size of 0.05:
<br />

```
MyESL.exe sample_files/angiosperm_100_sample_alns.txt --classes sample_files/test_pred.txt --DrPhylo --lambda1_grid 0.05,0.1,0.05 --lambda2_grid 0.05,0.1,0.05 --output sample_grid2x2_output 

```
<br />


### Model grid from DrPhylo

DrPhylo also generates outputs an ensemble clade model when the grid search option is used:

```

M-Grid_{clade_ID}.txt : Dataframe containing GSC, SCP, and CP.

M-Grid_{clade_ID}.png : Grid representation of the ensemble clade model using *_GSC_median.txt. This visualization will also contain SCP for all species in the clade of interest.
                                       The taxa with the lowest SCP will be at the top of the grid, and the lowest SCP for the clade is defined as the CP for the clade of interest. 

```
Note: The `{clade_ID}` will be replaced by the text file name when the `--classes` option is used for providing a user-defined phylogenetic hypothesis. In addition, DrPhylo generates summary files of Gene (GSS), Position (PSS), and Hypothesis (HSS) sparsity scores along with species prediction and probability score (SPS_SPP).

## Installation of MyESL into Python for DrPhylo analysis ##

To run DrPhylo, you will need Python 3.8 or later installed, as well as the following Python libraries:
```R
numpy
biopython
matplotlib
pandas
```

You can install these libraries using pip:

`pip install biopython numpy pandas matplotlib`

To perform DrPhylo analysis, get the ESL package using the following on the command line:

	git clone -b grid_search https://github.com/kumarlabgit/ESL ESL
	cd ESL
	bash setup.sh

## References ##
If you use DrPhylo in your research, please cite our articles:

1. Sharma, S. & Kumar, S. (2024). Discovering fragile clades and causal sequences in phylogenomics by evolutionary sparse learning. (In review)
2. Kumar, S. & Sharma, S (2021). Evolutionary Sparse Learning for Phylogenomics, Molecular Biology and Evolution, Volume 38, Issue 11, November 2021, Pages 4674–4682.
3. Sanderford, M., Sharma, S., Tamura, K., Stecher, G., Liu, J., Ye, J., and Kumar, S. (2024).  MyESL: A Software for Evolutionary Sparse Learning in Molecular Phylogenetics and Genomics (Submitted).
