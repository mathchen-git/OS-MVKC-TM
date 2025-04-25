% Multi-view Subspace Clustering on Topological Manifold
% function [result, Z, Tim, Wv, OBJ] = MVCsubspace_TM(data,labels, alpha, beta, knn, K_Lq, NITER)
function [result, Z, Tim, Wv, OBJ, H] = MVCsubspace_TM(data,labels, alpha, beta, knn, KP, Hq, NITER)
% data: cell array, view_num by 1, each array is num_samp by d_v
% num_clus: number of clusters
% num_view: number of views
% num_samp
% k: number of adaptive neighbours
% labels: groundtruth of the data, num by 1
%
% figure()
% [S1,~] = mapminmax(S, 0, 1)
% imshow(S1) 
% figure()
% [S1,~] = mapminmax(S, 0, 1)
% imshow(S1,'InitialMagnification','fit')
% colormap('jet');
%
if nargin < 3
    alpha = 1;
end
if nargin < 4
    beta = 1;
end
if nargin < 5
    knn = 20;
end
if nargin < 6
    normData = 1;
end
num_view = size(data,1);
num_samp = size(labels,1);
num_clus = length(unique(labels));
Wv = 1/num_view*ones(1,num_view);
% lambda = randperm(30,1);
lambda = 1;
% NITER = 30;
zr = 1e-10;

% =====================   Normalization =====================
% if normData == 1
%     for i = 1 :num_view
%         for  j = 1:num_samp
%             normItem = std(data{i}(j,:));
%             if (0 == normItem)
%                 normItem = eps;
%             end
%             data{i}(j,:) = (data{i}(j,:) - mean(data{i}(j,:)))/normItem;
%         end
%     end
% end
% NlzData = NormalizeData(normData,num_view,num_samp,data);

%  ====== Initialization =======
% claculate G_v for all the views
Zv = cell(num_view,1);
L = cell(num_view,1);
Lv = cell(num_view,1);
Dv = cell(num_view,1);
Zv2 = cell(num_view,1);
sumZ = zeros(num_samp);
for v = 1:num_view
%     Zv2 = constructW_PKN(NlzData{v}',knn);
%     Dv = diag(sum(Zv2));
%     Lv = Dv - Zv2;
%     Zv{v} = Zv2;
%     L{v} = Lv;
%     sumZ = sumZ + Zv2;
%     clear Zv Lv
%     Zv2{v} = constructW_PKN(NlzData{v}',knn);
    Zv2{v} = constructW_PKN(data{v}',knn);
    Dv{v} = diag(sum(Zv2{v}));
    Lv{v} = Dv{v} - Zv2{v};
    Zv{v} = Zv2{v};
    L{v} = Lv{v};
    sumZ = sumZ + Zv2{v};
end 
tic; 
% initialize S
% S = sumZ/num_view;

% initialize Z
Z = sumZ/num_view;


% initialize F
Z0 = Z-diag(diag(Z));
w0 = (Z0+Z0')/2;
D0 = diag(sum(w0));
L0 = D0 - w0;
[F0,~,~] = eig1(L0,num_clus,0);
F = F0(:,1:num_clus);

I = eye(num_samp);

% update ...
for iter = 1:NITER
    %
    % update Z_v
    [Zv,L_Zv] = updateZv(num_view,Wv,KP,Dv,Z,I);
    
    % update F
    [F,~,~] = updateF(Z,F,num_clus);
    
    % update Z
     Z = updateZ_K0_3(Z,num_samp,I,num_view,Wv,L,F,alpha,beta,KP,Hq);
    
    
    % update KP2  
    KP2 = updateK_Lq(num_samp,num_view,Z,Zv,alpha,Hq,I);
   
%       K_Lq = K_Lq./(Z+0.001); 
%     K_Lq = S'*K_Lq*S;



    
    % update w_v
    Wv = updateWv(num_view,Z,L_Zv);
   
     
    LZ = diag(sum(Z)) - Z;

sumobj_3=0;
for v = 1:num_view
        obj_3 = alpha*norm(Hq(v)*KP2{v}-Z,'fro')^2;
        sumobj_3=sumobj_3+obj_3;
end
obj_3_4 = sumobj_3+beta*trace(F'*LZ*F);

    for v = 1:num_view
        obj = obj_3_4 + trace(KP2{v}-2*KP2{v}*Zv{v}+Zv{v}'*KP2{v}*Zv{v}) + Wv(v)*trace(Z'*L_Zv{v}*Z);
    end
    OBJ(iter) = obj; 
    %
   

end
Tim = toc;
% =====================  result =====================
[groups, H] = SpectralClustering(Z,num_clus);
result = ClusteringMeasure(labels, groups);  % return metric

end


