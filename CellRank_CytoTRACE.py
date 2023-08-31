# CytoTrace approach

import anndata
import scanpy as sc
import scvelo as scv
import cellrank as cr
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys
import re


# In case the function refuses to write the file as a pdf. Reload image and convert to pdf
from PIL import Image


# Load in the Seurat-converted data (count-data) and the UMAP coordinates (as they are not always saved during this process)
hen_cd4 = sc.read_h5ad("/data/projects/Project_Gino/Trajectory_CD4_Amelie/Fixed_plots/sce.HEN.CD4_all_slinged_counts.h5ad")
X_UMAP = pd.read_csv("/data/projects/Project_Gino/Trajectory_CD4_Amelie/Additional_files/X_harmony_umap_coors.tsv", sep = "\t")


# The CytoTrace implementation in CellRank does not use the true spliced/unspliced count-table as the velocyto approach does, yet requires the keys in the scanpy object to work.
# Thus we create the required fields using the raw counts insead (according to the CellRank vignettes this was a method to to this)
df = pd.DataFrame(hen_cd4.X.toarray())
df
genes = hen_cd4.var['features'].to_list()
cells = hen_cd4.obs.index.to_list()

df.index = cells
df.columns = genes

hen_cd4.obs
hen_cd4.obsm['X_UMAP'] = X_UMAP.values


hen_cd4.layers["spliced"] = df[hen_cd4.var['features'].to_list()]
hen_cd4.layers["unspliced"] = df[hen_cd4.var['features'].to_list()]

# filter, normalize total counts and log-transform
sc.pp.filter_genes(hen_cd4, min_cells=10)
scv.pp.normalize_per_cell(hen_cd4)
sc.pp.log1p(hen_cd4)

df[hen_cd4.var['features'].to_list()]


# After the sc.pp.filter_genes a limit is set on the dimensions to be used
scv.pp.moments(hen_cd4, n_pcs=16, mode = 'distances', n_neighbors = 50, use_highly_variable = False, ) # might set enforce=True


# Use the CytoTrace kernel to determine pseudotime values for each cell
# As the initial pseudotime direction was opposite the the expected direction the pseudotime calculation was set to backwards
from cellrank.kernels import CytoTRACEKernel
ctk = CytoTRACEKernel(hen_cd4, backward=True)
ctk.compute_cytotrace(aggregation = 'median')
ctk.compute_transition_matrix(threshold_scheme="hard", nu=0.5, frac_to_keep=0.1)
ctk.plot_projection(basis="UMAP", color="ann_check_alt", legend_loc="right", dpi = 64, save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CytoTrace_plot_4.pdf")
image_1 = Image.open(r'/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CytoTrace_plot_4.png')
im_1 = image_1.convert('RGB')
im_1.save(r'/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CytoTrace_plot_4.pdf')


ctk.plot_random_walks(n_sims=1000,
                      start_ixs={"ann_check_alt" : "CD4_EM (Naive focus)"},
                      max_iter=0.35,
                      basis="UMAP",
                      color="ann_check_alt",
                      legend_loc="right",
                      seed=0,
                      n_jobs=100,
                      save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CytoTrace_plot_Random_Walk.pdf")


scv.pl.scatter(hen_cd4,
               c=["ct_pseudotime", "ann_check_alt"],
			   basis="UMAP",
			   legend_loc="right",
			   color_map="gnuplot2",
               color= 'ann_check_alt',
               palette= ["#F8766D", "#CF9400", "#E7861B", "#00B81F", "#00BC59"],
			   save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CytoTrace_plot_3.pdf"
			   )


# As it was set backwards it might have changed a bit
scv.pl.velocity_embedding_stream(hen_cd4,
                                 vkey= "T_bwd",
                                 basis = "UMAP",
                                 legend_loc="right",
                                 color= 'ann_check_alt',
                                 arrow_size=1,
                                 dpi = 64,
                                 palette= ["#F8766D", "#CF9400", "#E7861B", "#00B81F", "#00BC59"],
                                 title= "HEN CD4 Biopsies - scVelo: Velocity embedding (alt)",
                                 save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CytoTrace_plot_Embedding_alt_3.pdf")



#scv.pl.velocity_embedding_stream(hen_cd4, vkey="T_fwd",  legend_loc="right",  color= 'ann_check_alt', arrow_length=3, arrow_size=2, dpi = 128, title= "HEN CD4 Biopsies - scVelo: Velocity embedding (alt)", save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CytoTrace_plot_Embedding_alt_2.pdf")
sc.pl.violin(hen_cd4, 
		     keys=["ct_pseudotime"], 
			 groupby="ann_check_alt", 
			 rotation=90,
			 save = "/HEN_CD4_CytoTrace_plot_Violin.pdf"
			 )


# Create a combined kernel for downstream analysis: For this combined kernel half is based on the pseudotime, the other half is based on cell-cell similarity
# Both kernels are then combined to create a transition matrix
from cellrank.kernels import ConnectivityKernel

ck = ConnectivityKernel(hen_cd4).compute_transition_matrix()
ck

# Combine impact of both kernels to continue
combined_kernel = 0.5 * ctk + 0.5 * ck
combined_kernel
combined_kernel.compute_cytotrace(aggregation = 'median')
combined_kernel.compute_transition_matrix(threshold_scheme="hard", nu=0.5, frac_to_keep=0.1)
combined_kernel.plot_projection(basis="UMAP", color="ann_check_alt", legend_loc="right", dpi = 64, save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CombinedKernel_plot_4.pdf")


image_1 = Image.open(r'/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CombinedKernel_plot_4.png')
im_1 = image_1.convert('RGB')
im_1.save(r'/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CombinedKernel_plot_4.pdf')



combined_kernel.plot_random_walks(n_sims=100,
                                  start_ixs={"ann_check_alt" : "CD4_EM (Naive focus)"},
                                  max_iter=0.35,
                                  basis="UMAP",
                                  color="ann_check_alt",
                                  legend_loc="right",
                                  seed=0,
                                  n_jobs=100,
                                  save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/HEN_CD4_CombinedKernel_plot_Random_Walk.pdf")



# Steps below are to determine the Start and end points of the pseudotime lineage using the combined kernel
# 
# pk_new = cr.kernels.PseudotimeKernel.from_adata(hen_cd4, key="T_bwd")
g2 = GPCCA(combined_kernel)#, backward = True)
# Compute the Schur-decomposition and plot the real Schur-decomposition to see the number macrostates we can use 
g2.compute_schur(n_components=100)
g2.plot_spectrum(real_only=True,  figsize=(15, 10), save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g2_plot_Schur_Decomposition.pdf")
# NOTE: the real Schur-decomposition plot suggest a number of states, yet similar to an Elbow-plot this is not mandatory. Suggest to play around with this a bit.
g2.compute_macrostates(n_states=20, cluster_key="ann_check_alt")
g2.plot_macrostates(which="all", legend_loc="right", s=100, save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g_plot_macro_test_advanced_1.pdf")
g2.plot_macrostate_composition(key="ann_check_alt", figsize=(7, 4), save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g_plot_macro_test_advanced_2.pdf")
g2.plot_coarse_T(annotate=False, save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g_plot_macro_test_advanced_3.pdf")
g2.predict_terminal_states()
g2.plot_macrostates(which="terminal", legend_loc="right", s=100, save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g_plot_macro_test_advanced_4.pdf")
g2.set_terminal_states(states=["CD4_Th1_1", "CD4_Th17", "CD4_FH"])
g2.plot_macrostates(which="terminal", legend_loc="right", s=100, save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g_plot_macro_test_advanced_5.pdf")
g2.predict_initial_states()
g2.plot_macrostates(which="initial", s=100, save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g_plot_macro_test_advanced_6.pdf")
g2.set_initial_states(states=["CD4_EM (Naive focus)"])
g2.plot_macrostates(which="initial", legend_loc="right", s=100, save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g_plot_macro_test_advanced_7.pdf")
g2.compute_fate_probabilities()
g2.plot_fate_probabilities(legend_loc="right", save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g_plot_macro_test_advanced_8.pdf")
cr.pl.circular_projection(hen_cd4, keys="ann_check_alt", legend_loc="right", save = "/data/projects/Project_Gino/Trajectory_CD4_Amelie/Plots/CellRank/g_plot_macro_test_advanced_9.pdf")

