

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import anndata as ad


# should be the name of image data in adata
tissue_section = "CytAssist_FFPE_Human_Breast_Cancer"

# file path where outs data is located
file_path = "/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_visium/"


# read in svg results
gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
# these were not in the data
gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]
len(gene_list)

gene_list



### Read in adata ###

# read data
adata_visium_full = sc.read_visium(file_path)
# make unique
adata_visium_full.var_names_make_unique()
# get mitochondrial gene expression info
adata_visium_full.var["mt"] = adata_visium_full.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_visium_full, qc_vars=["mt"], inplace=True)


# make spatial position str to integer
# https://discourse.scverse.org/t/data-fomr-new-spatial-transcriptomics-from-10x/1107/6
# adata_visium_full.obsm['spatial'] = adata_visium_full.obsm['spatial'].astype(int)

# # get new coordiantes
# visium_aligned_coords = np.load("/home/caleb/Desktop/improvedgenepred/data/breastcancer_xenium_sample1_rep1/aligned_visium_points_to_xenium_image.npy")
# # add the new coordinates
adata_visium_full.obsm['spatial'] = adata_visium_full.obsm['spatial'].astype(int)

# add the image
# adata_visium_full.uns['spatial'] = img

# need to add this for subsetting
adata_visium_full.obs.index = adata_visium_full.obs.index.astype(str)

# log normalize the data
sc.pp.log1p(adata_visium_full)


# subet gene list
# adata_visium = adata_visium[:, gene_list]

# plot the data
sc.pl.spatial(adata_visium_full, color="TUBB2B", spot_size=150, title="TUBB2B", cmap="viridis")
# sc.pl.spatial(adata_visium_full, color="TUBB3", spot_size=150, title="TUBB3", cmap="viridis")
# sc.pl.spatial(adata_visium_full, color="TUBB4A", spot_size=150, title="TUBB4A", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="TUBB2A", spot_size=150, title="TUBB2A", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="TUBB8", spot_size=150, title="TUBB8", cmap="viridis")
# sc.pl.spatial(adata_visium_full, color="TUBB8B", spot_size=150, title="TUBB8B", cmap="viridis")
# sc.pl.spatial(adata_visium_full, color="TUBB2BP1", spot_size=150, title="TUBB2BP1", cmap="viridis")
# sc.pl.spatial(adata_visium_full, color="TUBB7P", spot_size=150, title="TUBB7P", cmap="viridis")



sc.pl.spatial(adata_visium_full, color="CEACAM8", spot_size=150, title="CEACAM8", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="CEACAM7", spot_size=150, title="CEACAM7", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="CEACAM6", spot_size=150, title="CEACAM6", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="CEACAM1", spot_size=150, title="CEACAM1", cmap="viridis")


sc.pl.spatial(adata_visium_full, color="OR4K1", spot_size=150, title="OR4K1", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="CUX2", spot_size=150, title="CUX2", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="PSG5", spot_size=150, title="PSG5", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="PDE11A", spot_size=150, title="PDE11A", cmap="viridis")

sc.pl.spatial(adata_visium_full, color="RAB30", spot_size=150, title="RAB30", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="FBXW7", spot_size=150, title="FBXW7", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="MTSS1", spot_size=150, title="MTSS1", cmap="viridis")

sc.pl.spatial(adata_visium_full, color="CEACAM6", spot_size=150, title="CEACAM6", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="CEACAM3", spot_size=150, title="CEACAM3", cmap="viridis")
sc.pl.spatial(adata_visium_full, color="CEACAM5", spot_size=150, title="CEACAM5", cmap="viridis")


########################################################################################################################



# should be the name of image data in adata
tissue_section = "CytAssist_FFPE_Human_Breast_Cancer"

# file path where outs data is located
file_path = "/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_visium/"


# read in svg results
gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
# these were not in the data
gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]
len(gene_list)

### Read in adata ###

# read data
adata_visium = sc.read_visium(file_path)
# make unique
adata_visium.var_names_make_unique()
# get mitochondrial gene expression info
adata_visium.var["mt"] = adata_visium.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_visium, qc_vars=["mt"], inplace=True)


# make spatial position str to integer
# https://discourse.scverse.org/t/data-fomr-new-spatial-transcriptomics-from-10x/1107/6
# adata_visium.obsm['spatial'] = adata_visium.obsm['spatial'].astype(int)

# # get new coordiantes
# visium_aligned_coords = np.load("/home/caleb/Desktop/improvedgenepred/data/breastcancer_xenium_sample1_rep1/aligned_visium_points_to_xenium_image.npy")
# # add the new coordinates
adata_visium.obsm['spatial'] = adata_visium.obsm['spatial'].astype(int)

# add the image
# adata_visium.uns['spatial'] = img

# need to add this for subsetting
adata_visium.obs.index = adata_visium.obs.index.astype(str)



import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc

# Step 1: Aggregate CEACAM6, CEACAM7, and CEACAM8
CEACAM_agg = adata_visium[:, ["CEACAM6", "CEACAM7", "CEACAM8", "CEACAM1"]].X.sum(axis=1)

# Step 2: Convert to sparse if necessary
if scipy.sparse.issparse(adata_visium.X):
    CEACAM_agg_sparse = scipy.sparse.csr_matrix(CEACAM_agg.reshape(-1, 1))
else:
    CEACAM_agg_sparse = CEACAM_agg.reshape(-1, 1)

# Step 3: Create new expression matrix by combining original X with CEACAM_agg
new_X = scipy.sparse.hstack([adata_visium.X, CEACAM_agg_sparse])

# Step 4: Update var to include the new gene
new_var = pd.concat(
    [adata_visium.var, pd.DataFrame(index=["CEACAM_agg"])],
    axis=0
)

# Step 5: Create the new AnnData object
adata_with_agg = sc.AnnData(X=new_X, var=new_var, obs=adata_visium.obs, obsm=adata_visium.obsm, uns=adata_visium.uns)

# Step 6: Verify the update
print(f"New AnnData shape: {adata_with_agg.shape}")  # Should reflect one extra gene
print(adata_with_agg.var.tail())  # Should show 'CEACAM_agg'



# log normalize the data
sc.pp.log1p(adata_with_agg)

# plot the data
sc.pl.spatial(adata_with_agg, color="CEACAM_agg", spot_size=150, title="CEACAM_agg", cmap="viridis")




########################################################################################################################



# should be the name of image data in adata
tissue_section = "CytAssist_FFPE_Human_Breast_Cancer"

# file path where outs data is located
file_path = "/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_visium/"


# read in svg results
gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
# these were not in the data
gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]
len(gene_list)

### Read in adata ###

# read data
adata_visium = sc.read_visium(file_path)
# make unique
adata_visium.var_names_make_unique()
# get mitochondrial gene expression info
adata_visium.var["mt"] = adata_visium.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_visium, qc_vars=["mt"], inplace=True)


# make spatial position str to integer
# https://discourse.scverse.org/t/data-fomr-new-spatial-transcriptomics-from-10x/1107/6
# adata_visium.obsm['spatial'] = adata_visium.obsm['spatial'].astype(int)

# # get new coordiantes
# visium_aligned_coords = np.load("/home/caleb/Desktop/improvedgenepred/data/breastcancer_xenium_sample1_rep1/aligned_visium_points_to_xenium_image.npy")
# # add the new coordinates
adata_visium.obsm['spatial'] = adata_visium.obsm['spatial'].astype(int)

# add the image
# adata_visium.uns['spatial'] = img

# need to add this for subsetting
adata_visium.obs.index = adata_visium.obs.index.astype(str)



import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc

# Step 1: Aggregate CEACAM6, CEACAM7, and CEACAM8
CEACAM_agg = adata_visium[:, ["RAB30", "FBXW7", "MTSS1"]].X.sum(axis=1)

# Step 2: Convert to sparse if necessary
if scipy.sparse.issparse(adata_visium.X):
    CEACAM_agg_sparse = scipy.sparse.csr_matrix(CEACAM_agg.reshape(-1, 1))
else:
    CEACAM_agg_sparse = CEACAM_agg.reshape(-1, 1)

# Step 3: Create new expression matrix by combining original X with CEACAM_agg
new_X = scipy.sparse.hstack([adata_visium.X, CEACAM_agg_sparse])

# Step 4: Update var to include the new gene
new_var = pd.concat(
    [adata_visium.var, pd.DataFrame(index=["CEACAM_agg"])],
    axis=0
)

# Step 5: Create the new AnnData object
adata_with_agg = sc.AnnData(X=new_X, var=new_var, obs=adata_visium.obs, obsm=adata_visium.obsm, uns=adata_visium.uns)

# Step 6: Verify the update
print(f"New AnnData shape: {adata_with_agg.shape}")  # Should reflect one extra gene
print(adata_with_agg.var.tail())  # Should show 'CEACAM_agg'



# log normalize the data
sc.pp.log1p(adata_with_agg)

# plot the data
sc.pl.spatial(adata_with_agg, color="CEACAM_agg", spot_size=150, title="gene_agg", cmap="viridis")





########################################################################################################################


### read in xenium data ###


# combined data
adata_xenium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/xeniumdata_visiumimage_data.h5ad')
adata_visium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/visiumdata_visiumimage_data.h5ad')

# make .X a csr matrix
adata_xenium.X = scipy.sparse.csr_matrix(adata_xenium.X)
adata_visium.X = scipy.sparse.csr_matrix(adata_visium.X)

# add array for gene expression
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)
adata_visium.X_array = pd.DataFrame(adata_visium.X.toarray(), index=adata_visium.obs.index)


# # log transform the data
# sc.pp.log1p(adata_xenium)
# sc.pp.log1p(adata_visium)


# plot total expression of genes
# adata_visium.obs["total_counts"] = adata_visium.X.sum(axis=1)
# adata_xenium.obs["total_counts"] = adata_xenium.X.sum(axis=1)

adata_visium.var['total_counts1'] = list(adata_visium.X.sum(axis=0).A1)
adata_xenium.var['total_counts1'] = list(adata_xenium.X.sum(axis=0).A1)

# plot the total counts of visium vs xenium as a scatter plot
plt.figure(figsize=(10, 5))
plt.scatter(np.log(adata_visium.var['total_counts1']), np.log(adata_xenium.var['total_counts1']))
plt.xlabel("Visium")
plt.ylabel("Xenium")
plt.title("Total counts of genes in visium vs xenium")
plt.xlim([min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5])
plt.ylim([min(np.log(adata_xenium.var['total_counts1']))-.5, max(np.log(adata_xenium.var['total_counts1']))+.5])
# plot a line
plt.plot([min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5], [min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5], color='red')
plt.show()



# plot the total counts of visium vs xenium as a scatter plot with gene names
plt.figure(figsize=(10, 5))
plt.scatter(np.log(adata_visium.var['total_counts1']), np.log(adata_xenium.var['total_counts1']))

# add gene names
for i, gene in enumerate(adata_visium.var_names):
    plt.text(np.log(adata_visium.var['total_counts1'][i]), np.log(adata_xenium.var['total_counts1'][i]), gene, fontsize=8)

plt.xlabel("Visium")
plt.ylabel("Xenium")
plt.title("Total counts of genes in visium vs xenium")
plt.xlim([min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5])
plt.ylim([min(np.log(adata_xenium.var['total_counts1']))-.5, max(np.log(adata_xenium.var['total_counts1']))+.5])
plt.show()







# read in text file
gene_summary = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/hisat2/gene_summary_transcriptome.txt", sep="\t")
# gene_summary = pd.read_csv("/home/caleb/Desktop/improvedgenepred/probes/bowtie/gene_summary_transcriptome.txt", sep="\t", header=None)


# Create the dictionary
gene_dict = {}

for i in range(0, len(gene_summary['Mismatches'])):
    gene_summary['Mismatches'][i]
    if pd.isna(gene_summary['Mismatches'][i]):
        gene_dict[gene_summary['Gene'][i]] = []
    else:
        gene_list = list(set(gene_summary['Mismatches'][i].split(",")))  # Split by ',' and get unique genes
        gene_dict[gene_summary['Gene'][i]] = gene_list

# Print 
print(gene_dict)

# get gene list
# read data
adata_visium_full = sc.read_visium("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_visium/")
gene_list = list(adata_visium_full.var_names)

gene_dict_filtered = {}

for i in range(len(gene_summary['Mismatches'])):
    gene_name = gene_summary['Gene'][i]  # Key
    mismatches = gene_summary['Mismatches'][i]  # Value

    if pd.isna(mismatches):
        gene_dict_filtered[gene_name] = []
    else:
        # Split, get unique genes, and filter by `gene_list`
        gene_set = set(mismatches.split(","))
        filtered_genes = [gene for gene in gene_set if gene in gene_list]

        gene_dict_filtered[gene_name] = filtered_genes  # Store filtered mismatches

# Print or use the dictionary
print(gene_dict_filtered)

# make dict a df
gene_summary_filtered = pd.DataFrame(gene_dict_filtered.items(), columns=["Gene", "Mismatches"])
# get rid or rows with empty lists in Misamatches column
gene_summary_filtered = gene_summary_filtered[gene_summary_filtered["Mismatches"].apply(lambda x: len(x) > 0)]
gene_summary_filtered


# sort genes
gene_summary = gene_summary.sort_values('Gene', ascending=True)
# get difference between columns 1 and 2
gene_summary["diff"] = gene_summary["Total_Count"] - gene_summary["Matched_Count"]
# set gene as index
gene_summary = gene_summary.set_index('Gene')
# merge df
merged_df = adata_visium.var.join(gene_summary, how="left")

# make adata.var the merged
adata_visium.var = merged_df









# plot the scatterplot but color by diff
plt.figure(figsize=(10, 5))
plt.scatter(np.log(adata_visium.var['total_counts1']), np.log(adata_xenium.var['total_counts1']), c=np.log(adata_visium.var["diff"]+1))
plt.xlabel("Visium")
plt.ylabel("Xenium")
plt.title("Total counts of genes in visium vs xenium")
plt.xlim([min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5])
plt.ylim([min(np.log(adata_xenium.var['total_counts1']))-.5, max(np.log(adata_xenium.var['total_counts1']))+.5])
# plot a line
plt.plot([min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5], [min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5], color='red')
plt.colorbar()
plt.show()





# plot the scatterplot but color by binary diff (0 or not 0)
colors = ['gray' if pd.isna(diff) else ('red' if diff != 0 else 'blue') 
          for diff in adata_visium.var["diff"]]
plt.figure(figsize=(10, 5))
plt.scatter(np.log(adata_visium.var['total_counts1']), np.log(adata_xenium.var['total_counts1']), c=colors)
plt.xlabel("Visium")
plt.ylabel("Xenium")
plt.title("Total counts of genes in visium vs xenium")
plt.xlim([min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5])
plt.ylim([min(np.log(adata_xenium.var['total_counts1']))-.5, max(np.log(adata_xenium.var['total_counts1']))+.5])
# plot a line
plt.plot([min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5], [min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5], color='red')
plt.show()



# plot the scatterplot but color by total alignments
plt.figure(figsize=(10, 5))
plt.scatter(np.log(adata_visium.var['total_counts1']), np.log(adata_xenium.var['total_counts1']), c=np.log(adata_visium.var['Total_Count']))
plt.xlabel("Visium")
plt.ylabel("Xenium")
plt.title("Total counts of genes in visium vs xenium")
plt.xlim([min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5])
plt.ylim([min(np.log(adata_xenium.var['total_counts1']))-.5, max(np.log(adata_xenium.var['total_counts1']))+.5])
# plot a line
plt.plot([min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5], [min(np.log(adata_visium.var['total_counts1']))-.5, max(np.log(adata_visium.var['total_counts1']))+.5], color='red')
plt.colorbar()
plt.show()
