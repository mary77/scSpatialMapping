%% input files
geometry=dlmread('../ALL_files/geometry.txt', ' ', 1, 0);
data1=csvread('../ALL_files/binarized_bdtnp.csv', 1, 0);
data1 = data1';
data2=csvread('../ALL_files/dge_binarized_distMap.csv', 1, 1);
ori_data1 = data1;
ori_data2 = data2;
[nGenes, nBins] = size(data1);

celldist = squareform(pdist(geometry));
% prob used to help select unique realLoc in training data (using all 84
% genes). Since geometry of in situ gene is precise, this weighting scheme
% will be very conservative - only for breaking ties
thres=min(min(celldist+(eye(nBins)*1000)));

prob = exp(-celldist/thres);

cc = corr(data1, data2);

[Y,realloc] = max(prob*cc);

% celldist2 used for selecting genes in RNAseq data
celldist2 = celldist(realloc, realloc);

%%
kfold = 10;
% fold = mod(randperm(n2), kfold);
% fold = fold + 1;
% dlmwrite('10fold.txt', fold');

fold = dlmread('fold_real_split.txt');
%%

targetN = 20;
remove = zeros(nGenes-targetN, kfold);
k = 10;
for f = 1:10
    disp(['fold ', int2str(f)]);
    tr_ins = (fold ~= f);
    data2 = ori_data2;
    data2_tr = data2(:, tr_ins);
    celldist2_tr = celldist2(tr_ins, tr_ins);
    removed = ones(nGenes, 1) < 1;    
    for i = 1:length(remove)
        input = find(~removed);
        x = backwardGeneSelectionNoInsituCV(data2_tr, k, celldist2_tr);
        data2(x, :) = [];
        data2_tr(x, :) = [];
        remove(i, f) = input(x);
        removed(remove(i, f)) = true;
    end
end
%% genes that have been eliminated are stored in the order of removal
save removed_genes remove