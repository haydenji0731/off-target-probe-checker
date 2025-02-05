############################################################################################################

# Plot two images side by side with gene expression colored patches


############################################################################################################

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import pandas as pd
import scanpy as sc
from PIL import Image
from utils import rasterizeGeneExpression_topatches, combine_adata_patches, plotRaster, prepare_data
import pickle
import scipy


# # file name
# file_name = "breastcancer_xenium_sample1_rep2"
# # resolution
# resolution = 250
# # read in the data
# adata = sc.read_10x_h5('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep2/cell_feature_matrix.h5')

# # Load the full-resolution spatial data
# # cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep2/breastcancer_xenium_sample1_rep2_visium_high_res_STalign.csv.gz", index_col=0)
# # cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/xenium_cell_centroids_visium_high_res.csv")
# cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep2/breastcancer_xenium_sample1_rep2_fullresolution_STalign.csv.gz", index_col=0)
# cell_centers

# # Load the full-resolution image
# Image.MAX_IMAGE_PIXELS = None
# img_name = "Xenium_FFPE_Human_Breast_Cancer_Rep2_he_image"
# img = np.array(Image.open("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/" + file_name + "/" + img_name + ".tif"))
# # img = np.load("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/visium_high_res_image.npy")
# plt.imshow(img)

# # add .obs
# adata.obs = cell_centers
# # add .obsm
# adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy().astype(int)
# # add image
# adata.uns['spatial'] = img
# # need to add this for subsetting
# adata.obs.index = adata.obs.index.astype(str)


# # get rid of genes that aren't in visium
# gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
# gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
# gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]
# # subset the data
# adata = adata[:, gene_list]

# # make an array of the gene expression data
# adata.X_array = pd.DataFrame(adata.X.toarray(), index=adata.obs.index)

# # need to subset bc there are negative values
# adata = adata[adata.obs["y_centroid"] > 0]
# adata = adata[adata.obs["x_centroid"] > 0]

# # Extract patches using `extract_patches_from_centers` function
# adata_patches = rasterizeGeneExpression_topatches(img, adata, patch_size=resolution, aggregation='sum')
# len(adata_patches)

# # scale the data
# scaling_factor = 1
# for i in adata_patches:
#     # aligned_visium_dictionary[i].X_array = sc.pp.log1p(aligned_visium_dictionary[i].X_array * scaling_factor)
#     adata_patches[i].X = sc.pp.log1p(np.round(adata_patches[i].X * scaling_factor))

# # combine the adata patches
# combined_adata = combine_adata_patches(adata_patches, img)

# # Example call to plot patches based on a specific obs column
# plotRaster(img, adata_patches, color_by='total_expression')
# plotRaster(img, adata_patches, color_by='gene_expression', gene_name='LPL')



############################################################################################################


# read in svg results
gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
# these were not in the data
gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]



### read in aligned data ###

# combined data
adata_xenium = sc.read_h5ad('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_sample1_rep1_aligned_toxeniumimage/combined_aligned_xenium_raw.h5ad')
adata_visium = sc.read_h5ad('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/visium_data_full.h5ad')

# make .X a csr matrix
adata_xenium.X = scipy.sparse.csr_matrix(adata_xenium.X)
adata_visium.X = scipy.sparse.csr_matrix(adata_visium.X)

# add array for gene expression
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)
adata_visium.X_array = pd.DataFrame(adata_visium.X.toarray(), index=adata_visium.obs.index)

# patches
with open('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_sample1_rep1_aligned_toxeniumimage/aligned_visium_dictionary_raw.pkl', 'rb') as f:
    aligned_visium_dictionary = pickle.load(f)

with open('/home/caleb/Desktop/improvedgenepred/data/breastcancer_sample1_rep1_aligned_tovisiumimage/visium_patches_full.pkl', 'rb') as f:
    aligned_xenium_dictionary = pickle.load(f)

# plt.imshow(adata_visium.uns['spatial'])

# scale the data
scaling_factor = 1
for i in aligned_visium_dictionary:
    # aligned_visium_dictionary[i].X_array = sc.pp.log1p(aligned_visium_dictionary[i].X_array * scaling_factor)
    aligned_visium_dictionary[i].X = sc.pp.log1p(np.round(aligned_visium_dictionary[i].X * scaling_factor))
    # aligned_visium_dictionary[i].X = sc.pp.scale(np.round(aligned_visium_dictionary[i].X * scaling_factor))

# adata_visium.X_array = adata_visium.X_array * scaling_factor

# scale the data
scaling_factor = 1
for i in aligned_xenium_dictionary:
    # aligned_visium_dictionary[i].X_array = sc.pp.log1p(aligned_visium_dictionary[i].X_array * scaling_factor)
    aligned_xenium_dictionary[i].X = sc.pp.log1p(np.round(aligned_xenium_dictionary[i].X * scaling_factor))
    # aligned_xenium_dictionary[i].X = sc.pp.scale(np.round(aligned_xenium_dictionary[i].X * scaling_factor))


# log transform the data
sc.pp.log1p(adata_xenium)
sc.pp.log1p(adata_visium)



############################################################################################################


import matplotlib.patches as patches
import matplotlib.pyplot

def plot_xenium_with_centers(adata, gene_list, g, patch_centers, gene_expression = True, patch_size=100, if_vis=True):
    """
    Plots the expression of a specified gene over an image using patch centers.

    :param adata: The AnnData object containing expression data and the image in `uns['spatial']`.
    :param gene_list: List of genes available in the expression data.
    :param g: The specific gene to visualize.
    :param patch_centers: Array of shape (N, 2) indicating the centers of each patch.
    :param patch_size: Size of the square patch (length of one side).
    :param if_pred: Boolean indicating whether to use "Predicted" or "True" expression title.
    """
    # Ensure the image is available in adata's uns
    if 'spatial' not in adata.uns:
        raise ValueError("The image data is missing in `adata.uns['spatial']`.")

    # Ensure the specified gene is available
    if g not in gene_list:
        raise ValueError(f"Gene '{g}' not found in the provided gene list.")

    # Extract the image from adata's uns
    image = adata.uns['spatial']

    # Create a matplotlib figure with the image in the background
    fig, ax = plt.subplots()
    ax.imshow(image)

    # make array
    adata.X_array = pd.DataFrame(adata.X.toarray(), index=adata.obs.index)

    if gene_expression:
        # Normalize the expression data of the specified gene
        gene_idx = gene_list.index(g)
        values = np.array(adata.X_array.iloc[:, gene_idx])
        norm = plt.Normalize(values.min(), values.max())
    else:
        # Normalize the expression data of the specified gene
        gene_idx = gene_list.index(g)
        values = np.array(adata.diff[:, gene_idx])
        norm = plt.Normalize(values.min(), values.max())

    # Create a colormap for the gene expression data
    scalar_map = plt.cm.ScalarMappable(norm=norm)

    # Calculate the top-left corner for each patch
    half_patch_size = patch_size // 2
    top_left_corners = patch_centers - np.array([half_patch_size, half_patch_size])

    # Plot each rectangular patch with the gene expression value, skipping missing data
    for i, (top_left_x, top_left_y) in enumerate(top_left_corners):
        if np.isnan(values[i]):
            continue  # Skip patches without valid gene expression values

        # Map the expression value to a color
        square_color = scalar_map.to_rgba(values[i])

        # Create and add the rectangle patch
        # ax.add_patch(patches.Rectangle(
        #     (top_left_x, top_left_y), patch_size, patch_size,
        #     linewidth=1, edgecolor='none', facecolor=square_color, alpha=1
        # ))

        # ax.add_patch(patches.Circle(
        #     (top_left_x, top_left_y), patch_size,
        #     linewidth=1, edgecolor='none', facecolor=square_color, alpha=1
        # ))

        # Create and add the hexagon patch
        hexagon = patches.RegularPolygon(
            (top_left_x, top_left_y), numVertices=6, radius=patch_size / 2,
            orientation=np.pi / 6, linewidth=1, edgecolor='none', facecolor=square_color, alpha=1
        )
        ax.add_patch(hexagon)

    # Remove axis ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # Display the color bar
    scalar_map.set_array([])
    fig.colorbar(scalar_map, ax=ax, orientation='vertical')

    # Determine and set the plot title
    if if_vis:
        title_prefix = "Visium"
    else:
        title_prefix = "Xenium"
    ax.set_title(f"{title_prefix} Expression of {g}", fontsize=10)

    # Display the plot
    # plt.savefig(f"/home/caleb/Desktop/improvedgenepred/output.png")
    plt.show()



g = "PDZD2"

plot_xenium_with_centers(adata_visium, list(adata_visium.var.index), g,  adata_visium.obsm['spatial'], patch_size=300, if_vis=True)
plot_xenium_with_centers(adata_xenium, gene_list, g,  adata_xenium.obsm['spatial'], patch_size=300, if_vis=False)








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
# gene_summary = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/hisat2/gene_summary_transcriptome.txt", sep="\t")
gene_summary = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/bowtie/gene_summary_transcriptome.txt", sep="\t")
# gene_summary = pd.read_csv("/home/caleb/Desktop/off-target-probe-checker/hisat2_genome/gene_summary_genome.txt", sep="\t")

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
len(gene_summary_filtered)


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
