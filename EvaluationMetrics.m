% function res =  EvaluationMetrics(Y, predY)
% % res = [clustering accuracy (ACC), normalized mutual information (NMI), Purity, 
% % Precision, Recall, F-score, and adjusted rand index (ARI)]
% % Ref.: Multiview Consensus Graph Clustering, Kun Zhan, Feiping Nie, Jing Wang, and Yi Yang.
% 
% [acc, nmi, Pu] = ClusteringMeasure(Y, predY);
% 
% ARI = RandIndex(Y, predY);
%  
% [Fscore, Precision, Recall] = compute_f(Y, predY);
% 
% res = [acc, nmi, Pu, Fscore, Precision, Recall, ARI];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  modify by Cuiling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Chen

function res =  EvaluationMetrics(Y, predY)

[ACC,NMI,Purity] = ClusteringMeasure(Y, predY);

res = [ACC,NMI,Purity];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
