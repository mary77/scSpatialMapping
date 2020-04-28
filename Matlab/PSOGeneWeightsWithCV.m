function [gBestWeights, mindist_train, mindist_test, weights] = PSOGeneWeightsWithCV(data1, data2, k, celldist, realloc, train_test_split, nAgents, nIter, weights)


data2_train = data2(:, train_test_split);
data2_test = data2(:, ~train_test_split);
realloc_train = realloc(train_test_split);
realloc_test =  realloc(~train_test_split);
m1 = size(data1,1);
m2 = size(data2_train,1);

assert(m1 == m2, 'data1 and data2 must have the same number of rows');

if (nargin < 9)
    %weights = [rand(m1, nAgents-1) ones(m1, 1)];
    weights = rand(m1, nAgents);
    weights = weights * diag(1./sum(weights));
end

v = rand(m1, nAgents);
v = v * diag(1 ./sum(v));

pBestWeights = weights;

pBestDist = zeros(nAgents, 1);

for i = 1:nAgents
    w = sparse(1:m1, 1:m1, weights(:, i));
    pBestDist(i) = calcDist(w*data1, w*data2_train, realloc_train, celldist, k);
end

[gBestDist, I] = min(pBestDist);
gBestWeights = weights(:, I);
w = sparse(1:m1, 1:m1, gBestWeights);

dist_test = calcDist(w*data1, w*data2_test, realloc_test, celldist, k);
disp([0 gBestDist dist_test])

i = 0;
while (i < nIter)
    i = i + 1;
    
    v = .2 * (pBestWeights - weights) * diag(rand(nAgents, 1)) + .2 * (repmat(gBestWeights, 1, nAgents) - weights) * diag(rand(nAgents, 1));
    
    disp([min(v(:)) prctile(v(:), 25) prctile(v(:), 50) prctile(v(:), 75) max(v(:))])
    
    weights = weights + v;
    
    weights (weights < 0) = 0;
    
    weights = weights * diag(1./sum(weights));
    
    currentDist = zeros(nAgents, 1);
    
    for j = 1:nAgents
        w = sparse(1:m1, 1:m1, weights(:, j));
        currentDist(j) = calcDist(w*data1, w*data2_train, realloc_train, celldist, k);
    end
    
    TF = currentDist < pBestDist;
    
    if (any(TF))
        pBestDist(TF) = currentDist(TF);
        pBestWeights(:, TF) = weights(:, TF);
        [Y, I] = min(pBestDist);
        if (Y < gBestDist)
            gBestDist = Y;
            gBestWeights = pBestWeights(:, I);
            w = sparse(1:m1, 1:m1, gBestWeights);
            dist_test = calcDist(w*data1, w*data2_test, realloc_test, celldist, k);
            disp([i gBestDist dist_test])
        else
            disp('o');
        end
    else
        disp(i);
    end
end

w = sparse(1:m1, 1:m1, gBestWeights);
mindist_train = calcDist(w*data1, w*data2_train, realloc_train, celldist, k);
mindist_test = calcDist(w*data1, w*data2_test, realloc_test, celldist, k);
end


function dist = calcDist(data1, data2, realloc, celldist, k)
cc = 1-pdist2(data1', data2', 'correlation');
cc(isnan(cc)) = 0;
[Y, I2] = sort(-cc);
[Y, I2] = sort(I2);
neigh = I2 <= k;

[I, J] = find(neigh);
J = realloc(J)';
ind=sub2ind(size(celldist), I, J);
aa = celldist(ind);
dist=mean(mean(reshape(aa, k, [])));
end