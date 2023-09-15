# distogram
Python script for calculating ensemble weighted inter-residue distance histogram (distogram)

## General description

Appropriate representation and comparison of ensembles of structural models are a topic of ongoing discussions in the field of structural biology.<sup>1,2,3</sup> Many traditionally used
ways of comparing the ensembles rely on RMSD-based superimposition of structures, which is very limited in its usefullness, especially in the case of flexible, multi-domain or
intrinsically disordered proteins. Superimposition-free model representations and scores have been under developement.<sup>1,4</sup> Inter-residue distance histograms (distograms)
are shown to be particularly promising, since they are: (1) superimposition-free (2) they are a model representation that can be robustly recovered (i.e. it is prior-independent)
in the ensemble reweighting procedures such as Maximum Entropy Method (3) it captures the global architecture of a molecule, for example relative orientation of domains to one another.<sup>1</sup>
The latter is one of the main advantages over other structure comparison methods, such as LDDT and GTS scores, which capture local accuracy of the model/assessment unit, but not the global
arrangement of subunits.<sup>4</sup>


For each of the conformational models _j_ in the ensemble, distances between all pair-wise combinations of C<sub>&alpha;</sub> atoms are computed.
By using  population fractions (weights) of ensemble members, C<sub>&alpha;</sub> inter-residue distance matrices are transformed into ensemble-weighted distograms, where each
point in the distogram represents an ensemble occupancy of a given distance bin. 



![distogram_illustration](https://github.com/mpopara/distogram/assets/40856779/3a41ee97-b559-4fe8-82a2-2095a3a55b1f)

Illustration is adapted from Dittrich _et al_, 2023.<sup>1</sup>

Counts (i.e. distance bin occupancies) and distance bins are saved as two separate files as numpy objects with extension .npy. 
Since distograms are multidimensional arrays, we can only visualize them in some kind of condensed form- for instance, by extracting the features/moments of distrograms.
Therefore, in the next step, moments of distograms are computed (mean, standard deviation and skewness), and after normalization, they are stacked into a R,G,B array (0-1 float), 
in order to create a RGB image, which will be exported as png and svg file.


## Input file requirements

* ensemble of conformational models provided as trajectory in any of the mdtraj compatible formats (dcd, nc, xtc..)
* topology as .pdb file
* .dat file containing weights (population fractions) of ensemble members. This space-delimited file is of a size N<sub>conformers</sub> x 2, where the first column contains indices of the ensemble members,
 and the second column contains their corresponding weights. This script assumes that the order of ensemble members in the trajectory file follows the same order as in the weights file.

## Dependencies
_distogram.py_ is a python script built on Python 3.8.8. Script was tested under the following configuration:

* Windows 10
* Python 3.8.8
* mdtraj 1.9.4
* numpy 1.23.0
* scipy 1.10.0
* numba 0.53.1
* matplotlib 3.7.1


## References
1. Dittrich, J.; Popara, M.; Kubiak, J.; Dimura, M.; Schepers, B.; Verma, N.; Schmitz, B.; Dollinger, . P.; Kovacic, F.; Jaeger, K. E.;
Seidel, C. A. M.; Peulen, T. O.; Gohlke, H., Resolution of Maximum Entropy Method-Derived Posterior Conformational Ensembles of a Flexible System Probed by FRET and Molecular Dynamics Simulations.
J Chem Theory Comput 2023, 19 (8), 2389-2409.

2. Tiberti, M.; Papaleo, E.; Bengtsen, T.; Boomsma, W.; Lindorff-Larsen, K., ENCORE: Software for Quantitative Ensemble Comparison. PLoS Comput Biol 2015, 11 (10)

3. Kufareva, I.; Abagyan, R., Methods of protein structure comparison. Methods Mol Biol 2012, 857, 231-57.

4. Mariani, V., Biasini, M., Barbato, A. & Schwede, T. lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests. Bioinformatics 29, 2722–2728 (2013)



## Authors

* Milana Popara
* Thomas-Otavio Peulen

