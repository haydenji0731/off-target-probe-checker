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


# file name
file_name = "breastcancer_xenium_sample1_rep2"
# resolution
resolution = 250
# read in the data
adata = sc.read_10x_h5('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep2/cell_feature_matrix.h5')

# Load the full-resolution spatial data
# cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep2/breastcancer_xenium_sample1_rep2_visium_high_res_STalign.csv.gz", index_col=0)
# cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/xenium_cell_centroids_visium_high_res.csv")
cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep2/breastcancer_xenium_sample1_rep2_fullresolution_STalign.csv.gz", index_col=0)
cell_centers

# Load the full-resolution image
Image.MAX_IMAGE_PIXELS = None
img_name = "Xenium_FFPE_Human_Breast_Cancer_Rep2_he_image"
img = np.array(Image.open("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/" + file_name + "/" + img_name + ".tif"))
# img = np.load("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/visium_high_res_image.npy")
plt.imshow(img)

# add .obs
adata.obs = cell_centers
# add .obsm
adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy().astype(int)
# add image
adata.uns['spatial'] = img
# need to add this for subsetting
adata.obs.index = adata.obs.index.astype(str)


# get rid of genes that aren't in visium
gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]
# subset the data
adata = adata[:, gene_list]

# make an array of the gene expression data
adata.X_array = pd.DataFrame(adata.X.toarray(), index=adata.obs.index)

# need to subset bc there are negative values
adata = adata[adata.obs["y_centroid"] > 0]
adata = adata[adata.obs["x_centroid"] > 0]

# Extract patches using `extract_patches_from_centers` function
adata_patches = rasterizeGeneExpression_topatches(img, adata, patch_size=resolution, aggregation='sum')
len(adata_patches)

# scale the data
scaling_factor = 1
for i in adata_patches:
    # aligned_visium_dictionary[i].X_array = sc.pp.log1p(aligned_visium_dictionary[i].X_array * scaling_factor)
    adata_patches[i].X = sc.pp.log1p(np.round(adata_patches[i].X * scaling_factor))

# combine the adata patches
combined_adata = combine_adata_patches(adata_patches, img)

# Example call to plot patches based on a specific obs column
plotRaster(img, adata_patches, color_by='total_expression')
plotRaster(img, adata_patches, color_by='gene_expression', gene_name='LPL')



############################################################################################################


# read in svg results
gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
# these were not in the data
gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]


### read in aligned data ###

# combined data
adata_xenium = sc.read_h5ad('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_sample1_rep1_aligned_toxeniumimage/combined_aligned_xenium_raw.h5ad')
adata_visium = sc.read_h5ad('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_sample1_rep1_aligned_toxeniumimage/combined_aligned_visium_raw.h5ad')

# make .X a csr matrix
adata_xenium.X = scipy.sparse.csr_matrix(adata_xenium.X)
adata_visium.X = scipy.sparse.csr_matrix(adata_visium.X)

# add array for gene expression
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)
adata_visium.X_array = pd.DataFrame(adata_visium.X.toarray(), index=adata_visium.obs.index)

# patches
with open('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_sample1_rep1_aligned_toxeniumimage/aligned_visium_dictionary_raw.pkl', 'rb') as f:
    aligned_visium_dictionary = pickle.load(f)

with open('/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_sample1_rep1_aligned_toxeniumimage/aligned_xenium_dictionary_raw.pkl', 'rb') as f:
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

# choose method
# method = "visium"
method = "xenium"

# resolution based on image
resolution = 250

# prepare the datasets
if method == "visium":
    X_train_tensor, y_train_tensor, scaled_coords, correct_order = prepare_data(aligned_visium_dictionary)
else:
    X_train_tensor, y_train_tensor, scaled_coords, correct_order = prepare_data(aligned_xenium_dictionary)



############################################################################################################



def plotRaster(ax, image, adata_patches, patch_size = 300, color_by='gene_expression', gene_name=None):
    """
    Plots patches on the original image, colored by either gene expression or a column in adata_patches.obs.
    
    Parameters:
    - ax: Matplotlib axis to plot on.
    - image: The original image array.
    - adata_patches: Dictionary of AnnData objects representing the patches.
    - color_by: How to color the patches ('gene_expression' or 'total_expression').
    - gene_name: The name of the gene to use if color_by is 'gene_expression'.
    """
    # Check inputs
    if color_by == 'gene_expression' and gene_name is None:
        raise ValueError("You must specify a gene_name when color_by='gene_expression'.")

    # Collect all values for normalization
    values = []
    for adata_patch in adata_patches.values():
        if color_by == 'gene_expression':
            expression = adata_patch.X[:, adata_patch.var_names.get_loc(gene_name)].sum()
            values.append(expression)
        elif color_by == 'total_expression':
            total_expression = adata_patch.X.sum()
            values.append(total_expression)
    
    # Get min and max values
    values = np.array(values)
    min_value, max_value = values.min(), values.max()

    # Plot the original image
    ax.imshow(image)

    # Plot each patch with the appropriate color
    for adata_patch in adata_patches.values():
        x_start, x_end, y_start, y_end = adata_patch.uns['patch_coords']
        
        if color_by == 'gene_expression':
            expression = adata_patch.X[:, adata_patch.var_names.get_loc(gene_name)].sum()
            normalized_value = (expression - min_value) / (max_value - min_value)
            color = plt.cm.viridis(normalized_value)
        elif color_by == 'total_expression':
            total_expression = adata_patch.X.sum()
            normalized_value = (total_expression - min_value) / (max_value - min_value)
            color = plt.cm.viridis(normalized_value)
        
        # Draw a rectangle for the patch
        # rect = mpatches.Rectangle((x_start, y_start), x_end - x_start, y_end - y_start,
        #                           linewidth=1, edgecolor='none', facecolor=color, alpha=1)
        # ax.add_patch(rect)

        # Create and add the hexagon patch
        hexagon = patches.RegularPolygon(
            (x_start, y_start), numVertices=6, radius=patch_size / 2,
            orientation=np.pi / 6, linewidth=1, edgecolor='none', facecolor=color, alpha=1
        )
        ax.add_patch(hexagon)

    # Create a color bar
    norm = plt.Normalize(min_value, max_value)
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.03, pad=0.04)
    cbar.set_label(f'{gene_name} Expression' if color_by == 'gene_expression' else "total_expression")

    ax.axis('off')

# Function to create combined plot for each gene
def plot_and_save_combined(gene_name, adata_visium, adata_xenium, image_visium, image_xenium,  save_path):
    """
    Plots both the gene expression raster and the scatter plot for each gene and saves them as one combined image.
    
    Parameters:
    - gene_name: Name of the gene to plot.
    - adata_visium: Visium patch data.
    - adata_xenium: Xenium patch data.
    - visimg_diff: Visium image difference data.
    - xenimg_diff: Xenium image difference data.
    - image_visium: Original Visium image.
    - image_xenium: Original Xenium image.
    - save_path: Path to save the output image.
    """
    # Create figure and axes
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    
    # Plot Visium patch data
    plotRaster(ax1, image_visium, adata_visium, color_by='gene_expression', gene_name=gene_name)
    ax1.set_title(f'Xenium Rep 1 Gene Expression: {gene_name}')
    
    # Plot Visium patch data
    plotRaster(ax2, image_xenium, adata_xenium, color_by='gene_expression', gene_name=gene_name)
    ax2.set_title(f'Xenium Rep 2 Expression: {gene_name}')

    # Save the combined plot
    # plt.suptitle("Correlation difference (Xenium - Visium): " + str(round(xenium_corr.loc[xenium_corr["Gene"] == gene_name, "Pearson"].values[0] - visium_corr.loc[visium_corr["Gene"] == gene_name, "Pearson"].values[0], 4)))
    plt.savefig(save_path)
    plt.close()



for i in range(267, len(gene_list)):
# for i in range(258, len(corr_diff['Gene'])):
    g = gene_list[i]
    plot_and_save_combined(
        gene_name=g,
        adata_visium=aligned_xenium_dictionary,
        adata_xenium=adata_patches,
        image_visium=adata_xenium.uns['spatial'],
        image_xenium=adata_xenium.uns['spatial'],
        save_path=f"/home/caleb/Desktop/improvedgenepred/results/xenium_xenium2_gene_plots_xeniumimage/{i}_{g}.png"
    )



############################################################################################################


g = "LPL"

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


g = "CEACAM6"

# plot_xenium_with_centers(adata_visium, gene_list, g,  adata_visium.obsm['spatial'], patch_size=300, if_vis=True)
plot_xenium_with_centers(adata_xenium, gene_list, g,  adata_xenium.obsm['spatial'], patch_size=300, if_vis=False)


