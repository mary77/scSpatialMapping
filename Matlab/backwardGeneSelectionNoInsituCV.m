function [remove, mindist] = backwardGeneSelectionNoInsituCV(data2_train, k, celldist2_train)

[m, n] = size(data2_train);

keep = ones(m, 1) > 0;
mindist = inf;

remove = [];
for i = 1:m
    keep(i) = false;
 
    cc2 = corr(data2_train(keep, :));
    [Y, I2] = sort(-cc2);
    [Y, I2] = sort(I2);
    neigh = I2 <= k;
    cc2=celldist2_train(neigh);
    % internal consistency in RNAseq
    cc2 = reshape(cc2, k, []);

    dist4 = mean(mean(cc2));
    
    if (dist4 < mindist)
        mindist = dist4;
        remove = i;
    end
    keep(i) = true;
end
