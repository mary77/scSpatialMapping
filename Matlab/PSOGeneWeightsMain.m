%%
geometry=dlmread('../ALL_files/geometry.txt', ' ', 1, 0);
celldist = squareform(pdist(geometry));

data1=csvread('../ALL_files/binarized_bdtnp.csv', 1, 0);
data1 = data1';
data2=csvread('../ALL_files/dge_binarized_distMap.csv', 1, 1);

ori_data1 = data1;
ori_data2 = data2;
cc = 1-pdist2(data1', data2', 'correlation');
[Y,realloc] = max(cc);

% top-k locations used for evaluation
k = 10;
[~, nCells] = size(data2);
%% 10-fold cv: save the split for easier reproducibility

fold = mod(randperm(nCells), 10) + 1;

dlmwrite('10fold.txt', fold);

%%
% fold = dlmread('10fold.txt'); 

%%
load removed_genes.mat remove
%%  optimal weight for 60 genes

disp([int2str(60), ' genes']);

nRemove = 84 - 60;

best_weights_60 = zeros(60, 10);
nAgents = 200;

for f = 1:10
    disp(['fold ', int2str(f)]);
    tr_ins = (fold~=f);
    removed = zeros(1,84)>0;
    removed(remove(1:nRemove, f)) = true;
    best_weights_60(:, f) = PSOGeneWeightsWithCV(data1(~removed,:), data2(~removed,:), k, celldist, realloc, tr_ins, nAgents, 40);
end
%% optimal weight for 40 genes

disp([int2str(40), ' genes']);
nRemove = 84 - 40;
best_weights_40 = zeros(40, 10);
nAgents = 200;

for f = 1:10
    disp(['fold ', int2str(f)]);
    tr_ins = (fold~=f);
    removed = zeros(1,84)>0;
    removed(remove(1:nRemove, f)) = true;
    best_weights_40(:, f) = PSOGeneWeightsWithCV(data1(~removed,:), data2(~removed,:), k, celldist, realloc, tr_ins, nAgents, 40);
end
%% optimal weight for 20 genes

disp([int2str(20), ' genes']);
nRemove = 84 - 20;
best_weights_20 = zeros(20, 10);
nAgents = 200;

for f = 1:10
    disp(['fold ', int2str(f)]);
    tr_ins = (fold~=f);
    removed = zeros(1,84)>0;
    removed(remove(1:nRemove, f)) = true;
    best_weights_20(:, f) = PSOGeneWeightsWithCV(data1(~removed,:), data2(~removed,:), k, celldist, realloc, tr_ins, nAgents, 40);
end
%%
save pso_gene_weights.mat best_weights_20 best_weights_40 best_weights_60