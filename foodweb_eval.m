%%
% Reproduce Table S10 in the supplementary material

%%
% Read data
data = dlmread('foodweb-clusters.txt', ',', 1, 0);
motif_embed = data(:, 2);
motif_rec   = data(:, 3);
edge_embed  = data(:, 4);
edge_rec    = data(:, 5);
louvain     = data(:, 6);
infomap     = data(:, 7);
gt1         = data(:, 8);
gt2         = data(:, 9);

%%
% Clustering results

% Classification 1
[AR_motif_embed1, F1_motif_embed1, ~, ~, ~] = RandIndex(gt1, motif_embed);
[AR_motif_rec1,   F1_motif_rec1,   ~, ~, ~] = RandIndex(gt1, motif_rec);
[AR_edge_embed1,  F1_edge_embed1, ~, ~, ~]  = RandIndex(gt1, edge_embed);
[AR_edge_rec1,    F1_edge_rec1, ~, ~, ~]    = RandIndex(gt1, edge_rec);
[AR_louvain1,     F1_louvain1, ~, ~, ~]     = RandIndex(gt1, louvain);
[AR_infomap1,     F1_infomap1, ~, ~, ~]     = RandIndex(gt1, infomap);
NMI_motif_embed1 = NMI(gt1, motif_embed);
NMI_motif_rec1   = NMI(gt1, motif_rec);
NMI_edge_embed1  = NMI(gt1, edge_embed);
NMI_edge_rec1    = NMI(gt1, edge_rec);
NMI_louvain1     = NMI(gt1, louvain);
NMI_infomap1     = NMI(gt1, infomap);
Pur_motif_embed1 = Purity(gt1, motif_embed);
Pur_motif_rec1   = Purity(gt1, motif_rec);
Pur_edge_embed1  = Purity(gt1, edge_embed);
Pur_edge_rec1    = Purity(gt1, edge_rec);
Pur_louvain1     = Purity(gt1, louvain);
Pur_infomap1     = Purity(gt1, infomap);

% Classification 2
[AR_motif_embed2, F1_motif_embed2, ~, ~, ~] = RandIndex(gt2, motif_embed);
[AR_motif_rec2,   F1_motif_rec2,   ~, ~, ~] = RandIndex(gt2, motif_rec);
[AR_edge_embed2,  F1_edge_embed2, ~, ~, ~]  = RandIndex(gt2, edge_embed);
[AR_edge_rec2,    F1_edge_rec2, ~, ~, ~]    = RandIndex(gt2, edge_rec);
[AR_louvain2,     F1_louvain2, ~, ~, ~]     = RandIndex(gt2, louvain);
[AR_infomap2,     F1_infomap2, ~, ~, ~]     = RandIndex(gt2, infomap);
NMI_motif_embed2 = NMI(gt2, motif_embed);
NMI_motif_rec2   = NMI(gt2, motif_rec);
NMI_edge_embed2  = NMI(gt2, edge_embed);
NMI_edge_rec2    = NMI(gt2, edge_rec);
NMI_louvain2     = NMI(gt2, louvain);
NMI_infomap2     = NMI(gt2, infomap);
Pur_motif_embed2 = Purity(gt2, motif_embed);
Pur_motif_rec2   = Purity(gt2, motif_rec);
Pur_edge_embed2  = Purity(gt2, edge_embed);
Pur_edge_rec2    = Purity(gt2, edge_rec);
Pur_louvain2     = Purity(gt2, louvain);
Pur_infomap2     = Purity(gt2, infomap);

%%
% table
results1 = [...
  AR_motif_embed1,  AR_motif_rec1,  AR_edge_embed1,  AR_edge_rec1,  AR_louvain1,  AR_infomap1; ...
  F1_motif_embed1,  F1_motif_rec1,  F1_edge_embed1,  F1_edge_rec1,  F1_louvain1,  F1_infomap1; ...    
  NMI_motif_embed1, NMI_motif_rec1, NMI_edge_embed1, NMI_edge_rec1, NMI_louvain1, NMI_infomap1; ...  
  Pur_motif_embed1, Pur_motif_rec1, Pur_edge_embed1, Pur_edge_rec1, Pur_louvain1, Pur_infomap1;];

results2 = [...
  AR_motif_embed2,  AR_motif_rec2,  AR_edge_embed2,  AR_edge_rec2,  AR_louvain2,  AR_infomap1; ...
  F1_motif_embed2,  F1_motif_rec2,  F1_edge_embed2,  F1_edge_rec2,  F1_louvain2,  F1_infomap1; ...    
  NMI_motif_embed2, NMI_motif_rec2, NMI_edge_embed2, NMI_edge_rec2, NMI_louvain2, NMI_infomap1; ...  
  Pur_motif_embed2, Pur_motif_rec2, Pur_edge_embed2, Pur_edge_rec2, Pur_louvain2, Pur_infomap1;];

[results1; results2]
    