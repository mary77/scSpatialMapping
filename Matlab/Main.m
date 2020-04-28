
%% first, perform gene selection to get 60, 40, and 20 genes, using SUPERVISED feature selection
GeneSelectionMain

%% Use PSO to obtain optimal gene weights for different folds of training data
PSOGeneWeightsMain

%% Make predictions on testing data using weights learned from above
makePredictions_subc1

makePredictions_subc2

makePredictions_subc3


