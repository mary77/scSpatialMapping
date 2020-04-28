Maryam Zand & Dr. Jianhua Ruan
Nov 2018

Abstact:
Spatial mapping of single cells in the Drosophila embryo from transcriptomic data based on topological consistency

Script to run the program:
main.m

This path should be added to matlab path :
"../ALL_files"
All_files contains scRNAseq and insitu expression data that are accessible from this link: 
https://www.synapse.org/#!Synapse:syn15665609/wiki/582909

Functions

%% first, perform gene selection to get 60, 40, and 20 genes, using SUPERVISED feature selection
GeneSelectionMain

%% Use PSO to obtain optimal gene weights for different folds of training data
PSOGeneWeightsMain

%% Make predictions on testing data using weights learned from above
makePredictions_subc1

makePredictions_subc2

makePredictions_subc3