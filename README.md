# Agent-Based Model Described in Rouillard et al. Journal of Physiology 2012
Our group developed an agent-based model of myocardial scar healing, implemented in MATLAB. The model tracks migration, proliferation, and collagen remodeling by fibroblasts over 6 weeks of simulated infarct healing, incorporating the influence of regional strains, chemokine gradients, and local matrix orientation on the orientation of the fibroblasts and the collagen fibers they produce. The original model was published in the Journal of Physiology in 2012, and has been used with various minor modifications in several other papers and abstracts from our group since that time. If you use or adapt the model, please cite the original Journal of Physiology article: [Rouillard AD, Holmes JW. Mechanical regulation of fibroblast migration and collagen remodeling in healing myocardial infarcts, J Physiol, 590(18):4585-4602, Sep 2012](http://www.ncbi.nlm.nih.gov/pubmed/22495588).

This repository contains two .m file with the MATLAB code for the Journal of Physiology model, as well as several input workspaces and two sample output workspaces.

# ABM_20130514_JPhysiolModel.m
This is the agent-based model code used in Rouillard & Holmes (2012) "Mechanical regulation of fibroblast migration and collagen remodeling in healing myocardial infarcts." J Physiology 590: p4585-4602.
The code is currently setup to sequentially simulate healing of a circular infarct under uniaxial and biaxial strain environments as discussed in the paper. Options may be changed to simulate other infarct geometries, environments, or conditions.

## Input workspaces
* ABM_ConcField_Circle.mat - This is a Matlab file containing the chemokine concentration field required to run the ABM code for circular infarcts.
* ABM_ConcField_Ellipse2_Horizontal.mat - This is a Matlab file containing the chemokine concentration field required to run the ABM code for elliptical infarcts with elongation ratio = 2, oriented horizontally.
* ABM_ConcField_Ellipse2_Vertical.mat - This is a Matlab file containing the chemokine concentration field required to run the ABM code for elliptical infarcts with elongation ratio = 2, oriented vertically.
* ABM_ConcField_Ellipse5_Horizontal.mat - This is a Matlab file containing the chemokine concentration field required to run the ABM code for elliptical infarcts with elongation ratio = 5, oriented horizontally.
* ABM_ConcField_Ellipse5_Vertical.mat - This is a Matlab file containing the chemokine concentration field required to run the ABM code for elliptical infarcts with elongation ratio = 5, oriented vertically.

## Output workspaces
* ABM_20130514PAR_Biax_Circle_AD_R3_Ws0.16667_42days.mat - This is the Matlab results file for a simulation of circular infarct healing under biaxial strain.
* ABM_20130514PAR_UniaxC_Circle_AD_R3_Ws0.16667_42days.mat - This is the Matlab results file for a simulation of circular infarct healing under uniaxial (circumferentially-oriented) strain.

# ABM_20130514_ResultsSummary.m
This is a post-processing code used to generate four figures from the '.mat' results files listed above.
* Figure 1 contains spatial maps of collagen density (denoted by color) and orientation (denoted by quiver directions) after 3 weeks of simulated healing under biaxial or uniaxial strains.
* Figure 2 contains spatial maps of cell positions and orientations (denoted by quiver directions) after 3 weeks of simulated healing under biaxial or uniaxial strains.
* Figure 3 contains plots of collagen fiber (top row) and cell (bottow row) mean vector lengths, mean vector angles, and densities over the 42 day healing time course under both biaxial and uniaxial strain conditions.
* Figure 4 contains histograms of collagen fiber and cell orientations after 3 weeks of simulated healing under biaxial or uniaxial strains.

#### NOTE: Due to stochasticity of the model, repeat simulations will yield very slight variations in results. The model is setup to run 4 repeat simulations, and the "...ResultsSummary.m" file averages across those 4 simulations in the final figures. The original paper only included 1 simulation, so those plots (though qualitatively identical) were less smooth
