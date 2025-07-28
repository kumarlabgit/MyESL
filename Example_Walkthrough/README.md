# A walkthrough with MyESL #
In this example walkthrough, we explore key MyESL options using the provided Fungi dataset. Each step includes relevant MyESL commands along with a discussion of the corresponding outputs.

<details>
<summary><strong>Data Input Options</strong></summary>

### Model Building using sequence alignments and a phylogeny in newick format with clade ID 

```
MyESL.exe Fungi_data\aln.txt --tree Fungi_data\Fungi_T1_with_ID.nwk
```

This command will create a directory ``output`` that will contain all output files. The phylogeny has two clade IDs, ``Clade_X1`` and ``Control``. MyESL will produce two clade models for these two clades. In this clade, the default values for the site penalty ``(lambda1)`` and gene penalty ``(lambda2)`` will be used, which is ``0.1`` for both cases. The ``output`` directory contains model grids (M_grid) for the clades ``Clade_X1`` and ``Control``. 

![image](https://github.com/user-attachments/assets/62fa9fd9-dc92-41e1-a1c5-1a6c955c4082)

### A clade model for a specific clade 
We can build clade models for a specific clade using a list of clade IDs provided by a text file 

```
MyESL.exe Fungi_data\aln.txt --tree Fungi_data\Fungi_T1_with_ID.nwk  --clade_list Fungi_data\clade_Control.txt
```

### Building a clade model using a phylogenetic tree without clade IDs
Clade IDs can be generated for clades in the input phylogeny. In this case, the size of the clades (number of species in clades) can be defined by the users. 

```
MyESL.exe Fungi_data\aln.txt --tree Fungi_data\Fungi_T1_with_ID.nwk  --gen_clade_list 20,40
```

MyESL will assign clade IDs for clades that have a minimum of 20 species or a maximum of 40 species. 

### Building a clade model without phylogenetic tree input
You can provide a phylogenetic hypothesis for species groupings using a **tab-separated text file**, which serves as the response input for sparse learning.

Each line in the file should follow this format:

Where `<response_value>` can be:
- `+1` for species **inside** the clade of interest
- `-1` for species **outside** the clade
- `0` for species to be **excluded** from the analysis

This allows you to define clade membership manually, without relying on a Newick-format phylogenetic tree.

![image](https://github.com/user-attachments/assets/6c0be5c9-ce3c-4d4b-a84d-5bee0be15970)

```
MyESL.exe Fungi_data\aln.txt --classes Fungi_data\A_B_Hyp.txt 
```

</details>

<details>
<summary><strong>Model Options</strong></summary>
MyESL uses site and gene sparsity parameters equal to ``0.1`` by default. These parameters can be set by ``lambda1`` and ``lambda2``. 

```
MyESL.exe Fungi_data\aln.txt --classes Fungi_data\A_B_Hyp.txt --lambda1 0.1 --lambda2 0.2 

```

Selecting appropriate parameter values can be a challenging task. To address this, MyESL allows the use of a parameter grid, enabling the construction of multiple clade models across a range of values. The ``lambda1_grid`` and ``lambda2_grid`` options define the parameter ranges using a minimum value, a maximum value, and the step size for each grid.

```
MyESL.exe Fungi_data\aln.txt --classes Fungi_data\A_B_Hyp.txt --lambda1_grid 0.1,0.9,0.1 --lambda2_grid 0.1,0.9,0.1 

```
When using grid search options (`--lambda1_grid` and `--lambda2_grid`), MyESL builds multiple models across combinations of sparsity parameters. While this allows thorough exploration, it can be time-consuming, and not all models may be a good fit for the data.

To improve computational efficiency, MyESL supports **model skipping** using the following options:

- `--min_groups`: Skips models that include fewer than the specified number of groups.
- `--grid_rmse_cutoff`: Skips models with a root mean squared error (RMSE) or model fit score (MFS) above the given threshold.
- `--grid_acc_cutoff`: Skips models with accuracy below the specified cutoff.

These options help focus the analysis on well-fitting models and reduce unnecessary computation.

```
MyESL.exe Fungi_data\aln.txt --classes Fungi_data\A_B_Hyp.txt --lambda1_grid 0.1,0.9,0.1 --lambda2_grid 0.1,0.9,0.1 --min_groups 3
MyESL.exe Fungi_data\aln.txt --classes Fungi_data\A_B_Hyp.txt --lambda1_grid 0.1,0.9,0.1 --lambda2_grid 0.1,0.9,0.1 --grid_rmse_cutoff 0.5
MyESL.exe Fungi_data\aln.txt --classes Fungi_data\A_B_Hyp.txt --lambda1_grid 0.1,0.9,0.1 --lambda2_grid 0.1,0.9,0.1 --grid_acc_cutoff 0.95
```
These options can also be used to select the best set of sparsity parameter values. For example, we want to select models with those sparsity parameters that have a training accuracy of 95%. 

MyESL also offers to build a model with mono-level sparsity at the site level. This is useful when biological boundaries for datasets are not defined.

```
MyESL.exe Fungi_data\aln.txt --classes Fungi_data\A_B_Hyp.txt --lambda1 0.1 --no_group_penalty

```

MyESL performs logistic LASSO regression by default. However, MyESL also offers to perform sparse group LASSO regression with least squared loss when the response is a continuous value (e.g., body mass). 

```
MyESL.exe Fungi_data\aln.txt --classes Path\Continuous_response.txt --lambda1 0.1 --lambda2 0.2 --method leastr

```
</details>

<details>
<summary><strong>Class Balancing in MyESL</strong></summary>
Class imbalance can significantly affect the performance of supervised machine learning models. In the context of phylogenetic modeling, when one class—such as taxa **within** the focal clade—is underrepresented, the model may become biased toward the **majority class**. This often leads to **high overall accuracy** but poor **sensitivity and specificity** for the minority class, thereby reducing the model's ability to correctly classify taxa **inside or outside** the clade and to identify truly informative features (e.g., genes or sites). To address this issue, MyESL provides several class balancing strategies through the `--class_bal` directive, including **class weighting**, **upsampling**, **downsampling**, and novel **phylogeny-aware balancing**. The phylogenetic-aware class balancing is performed if the phylogenetic hypothesis is provided using a phylogeny with a clade ID. 

```
MyESL.exe Fungi_data\aln.txt --tree Fungi_data\Fungi_T1_with_ID.nwk  --class_bal phylo
```
Other class balancing can be performed for both the tree input and class input using a text file. 

```
MyESL.exe Fungi_data\aln.txt --tree Fungi_data\Fungi_T1_with_ID.nwk  --class_bal up
MyESL.exe Fungi_data\aln.txt --tree Fungi_data\Fungi_T1_with_ID.nwk  --class_bal down
MyESL.exe Fungi_data\aln.txt --tree Fungi_data\Fungi_T1_with_ID.nwk  --class_bal weighted
```

</details>

<details>
<summary><strong>DrPhylo with MyESL</strong></summary>
MyESL enables DrPhylo analysis for a specified clade, either defined by clade ID using a phylogeny in a NEWICK format or via a response file in text format. DrPhylo builds multiple sparse models across combinations of site and group sparsity parameters. By default, it explores values from 0.1 to 0.9 (in steps of 0.1), generating 81 models. Users can customize this range using the "--lambda_grid" directive. To reduce computation, DrPhylo skips models that include fewer than three genes. It outputs summary statistics, including PSS, GSS, and HSS scores, as well as a model grid. 

```
MyESL.exe Fungi_data\aln.txt --classes Path\Continuous_response.txt --lambda1 0.1 --lambda2 0.2 --DrPhylo

```

</details>

<details>
<summary><strong>MyESL Outputs</strong></summary>

Users can define an output directory name. This directory will be created in the current working directory. If no output directory name is defined by users, a directory named “output” will be created by default. 

```
MyESL.exe Fungi_data\aln.txt --classes Fungi_data\A_B_Hyp.txt --lambda1 0.1 --lambda2 0.2 --output Fungi_Clade_A_B

```
This will create an output directory of the name, ``Fungi_Clade_A_B``.

Users can customize the graphical outputs generated by MyESL. For example, in the **DrPhylo** analysis, the size of the **model grid output** can be adjusted. By default, the grid is set to **20×20** (20 rows for species and 20 columns for genes). To modify the grid size, use the `--m_grid <row, col>` directive.

```
MyESL.exe Fungi_data\aln.txt --classes Path\Continuous_response.txt --lambda1 0.1 --lambda2 0.2 --DrPhylo --m_grid 20,20

```

Users can also output different sparsity scores using the directive ``--stats_out``, where 
P: Site-level (position) sparsity scores
G: Gene or group sparsity scores
H: Hypothesis sparsity score
S: Species prediction scores

```
MyESL.exe Fungi_data\aln.txt --classes Path\Continuous_response.txt --lambda1 0.1 --lambda2 0.2 --DrPhylo --m_grid 20,20 --stats_out GS

```
This will output group sparsity scores and species prediction scores in text format. 

</details>

<details>
<summary><strong>Prediction using MyESL</strong></summary>
MyESL provides a separate pipeline `` MyESL_model_apply.exe `` for applying the ESL model to predict traits or determine the clade membership of new species. To use this feature, users must provide sequence alignments that are aligned with the training alignments. For example, we want to predict the membership of a list of species using the ESL model for the control clade. The key inputs are 

```
ESL model:      An ESL model built for the clade of interest. The ESL model is stored in a directory named starting with MyESL_model. 
Alignment list: A text file containing the full path for the alignments of new species. These alignments must be aligned with the sequence alignments of the test species. 
Response file:  A text file containing the name of the species (one per line) for which we want the prediction score. 
```
An example Model file 

An example response file 
<img width="278" height="191" alt="image" src="https://github.com/user-attachments/assets/7f9d0894-e78f-4308-917b-2c2dfde96ae0" />


</details>

  

