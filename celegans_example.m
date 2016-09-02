% Generate result from Figure 2 in
% Higher-order organization of complex networks. 
% Austin R. Benson, David F. Gleich, and Jure Leskovec.
% Science, 353.6295 (2016): 163--166.

load 'celegans_data';
W = MotifAdjacency(A, 'bifan');
[LCC, lcc_inds] = LargestConnectedComponent(W);
[cluster, condv, condc] = SpectralPartitioning(LCC);
orig_inds = lcc_inds(cluster);
neurons = labels(orig_inds);

% Plot in 2d space with labels
scatter(pos(:, 1), pos(:, 2), 30, 'black', 'filled'); hold on;
scatter(pos(orig_inds, 1), pos(orig_inds, 2), 30, 'red', 'filled');
for i = 1:length(orig_inds)
    text(pos(orig_inds(i), 1), ...
         pos(orig_inds(i), 2), ...
         neurons(i));
end
xlabel('x (mm)');
ylabel('y (mm)');
title('C. elegans frontal neurons bifan cluster');
