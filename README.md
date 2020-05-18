## Mouse Lifespan Synaptome Atlas code repository

Created by Zhen Qiu (z.qiu@ed.ac.uk), Erik Fransén (erikf@kth.se) and  Seth Grant (seth.grant@ed.ac.uk)


### Introduction
Matlab code to generate figures for the manuscript 'Cizeron, M.\*, Qiu, Z.\*, Koniaris, B., Gokhale, R., Komiyama, N., Fransén, E., Grant, S.G.N. (2020). A brain-wide atlas of synapses across the mouse lifespan. (* Equal contribution)'.  



### License
The code  is licensed  under the MIT License (refer to the LICENSE file for details).


### Prerequisites
1. Mathworks Matlab 2014b or above
2. Export_fig: A MATLAB toolbox for exporting publication quality figures https://github.com/altmany/export_fig
3. Matlab Toolbox for Bayesian Estimation (MBE) https://github.com/NilsWinter/matlab-bayesian-estimation


### Descriptions
#### 1. compare_3Mto18M_full_parameter
Matlab code to generate Figure. 1C: compare the puncta parameters of 3M with that of 18M on 109 brain regions.
#### 2. compare_3Mto18M_subtypes
Matlab code to generate Figure. S13: compare the puncta subtype densities of 37 subtype between 3M with that of 18M in 109 brain regions.
#### 3. compare_hypersimi_2Wto3M_2Wto18M 
Matlab code to test  if the 18M, in contrast to 3M, is more similart to 2W. (Figure S19, 18M2W versus 3M2W). This was done by using the hypersimilarity matrix Figure 3C.
#### 4. compare_hypersimi_2Wto3M_2Wto18M_HPF
Matlab code to test  if the 18M in Hippocampal formation, in contrast to 3M, is more similar to 2W. (Figure S19, 18M2W versus 3M2W). This was done by using the hypersimilarity matrix Figure 3D.
#### 5. permutation_test_simi_ratio, permutation_test_simi_ratio_time, permutation_test_simi_ratio_spacetime, 
Permutation test to test the significance of the similarity ratio. Firstly permute the regions of similarity matrices (Figure S16) in space for each age group, then permute the matrices (Figure S16) of different age groups of the same region for time, finally permute the hypersimilarity matrix (Figure 3D) in both space and time.
#### 6. plot_heatmap_trajectories_full_parameters
Matlab code to plot the heatmaps(Figure S5) of trajectories of PSD95 and SAP102 parameters. 
#### 7. plot_class_heatmap
Matlab code to plot the heatmaps (Figures 2A-C, S9-11) of 3 types and 37 subtypes on 12 brain regions.
#### 8. plot_diversity_lifespan
Matlab code to plot the diversity lifespan trajectories (Figure 2F) on 12 brain regions.
#### 9. plot_diversity_lifespan_unsupervised
Matlab code to plot the unsupervised diversity brain maps (Figures 2G, S15).
#### 10. plot_correlation_classpercent_simiratio
Matlab code to plot the correlation (Figure S20)between lifespan trajectories of synapse subtype percentage and similarity ratio on 12 brain regions .
#### 11. plot_hypersimi_matrix
Matlab code to plot the hypersimilarity matrix (Figures 3C and S18) for the whole section.
#### 12. plot_hypersimi_matrix_HPF
Matlab code to plot the hypersimilarity matrix (Figure 3D) of hippocampal formation .
#### 13. plot_ratio_withintobetweenregion
Matlab code to calculate and plot the lifespan trajectories of similarity ratio(Figure 3B).
#### 14. plot_ratio_withintobetweenregion_per_region
Matlab code to calculate and plot the lifespan trajectories of similarity ratio  in 12 brain regions (Figure S17).

#### 15. synaptomeModel
Matlab code to build the computational model (Figures 4B and S22) that illustrates how lifespan synaptome changes affect synaptic responses to distinct temporal patterns of neural activity 


### References

Cizeron, M.\*, Qiu, Z.\*, Koniaris, B., Gokhale, R., Komiyama, N., Fransén, E., Grant, S.G.N. (2020). A brain-wide atlas of synapses across the mouse lifespan. (* Equal contribution)
