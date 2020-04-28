%%
geometry=dlmread('../ALL_files/geometry.txt', ' ', 1, 0);
celldist = squareform(pdist(geometry));
nBins = size(geometry,1);
% sigma value for probability calculation
sigma=median(min(celldist+(eye(nBins)*1000)));

% converts euclidean distance to probabilities: pij = exp(-dij / sigma)
prob = exp(-celldist/sigma);

data1=csvread('../ALL_files/binarized_bdtnp.csv', 1, 0);
data1 = data1';
data2=csvread('../ALL_files/dge_binarized_distMap.csv', 1, 1);
nMarkerGenes = size(data1,1);
%% 40 genes to be used

load removed_genes remove
removed = zeros(1,nMarkerGenes)>0;
nRemove = nMarkerGenes - 40;
removed(remove(1:nRemove)) = true;
%%
load pso_gene_weights.mat best_weights_40
%%
fold = dlmread('10fold.txt');

% number of predicted locations per cell
k = 10;

predictions = zeros(k, 1297);
finalCC = zeros(3039,1297);
% make predictions for each fold using weights from corresponding training data
for f = 1:10
    disp(['fold ', int2str(f)]);
    test_ins = (fold==f);
    %w = diag(best_weights_40(:, f));
    %%d1 = w*data1(~removed, :);
    %d2 = w*data2(~removed, test_ins);
    d1 = data1(~removed, :);
    d2 = data2(~removed, test_ins);
    cc = 1-pdist2(d1', d2', 'correlation');
    % optional - compute weighted sum of predictions by neighbors
    cc = prob * cc;
    finalCC(:,test_ins) = cc;
    % top-10 locations for each cell
    [Y, I] = sort(-cc);
    [Y, I] = sort(I);
    neigh = I <= k;    
    [I, J] = find(neigh);
    top10predictions=reshape(I, 10, []);

    % sort the predictions
    aa = cc(neigh);
    aa = reshape(aa, k, []);
    [X, Y] = sort(-aa);
    top10predictions_ranked = [];
    for i = 1:size(top10predictions, 2)
        top10predictions_ranked = [top10predictions_ranked top10predictions(Y(:, i), i)];
    end
    
    predictions(:, test_ins) = top10predictions_ranked;    
end
finalCC = exp(finalCC);
%%
selected40Genes = find(~removed);
save subchallenge2_results.mat selected40Genes predictions


