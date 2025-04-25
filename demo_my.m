clear
clc

path = './';
addpath(genpath(path));
dataName = 'prokaryotic_uni'; 
% load([path,'dataset/',dataName],'X','Y');
% data = X';
% data = X;  % only for yaleA_3view

load([path,'dataset/',dataName],'XX','Y');
data = XX';

load([path,'ConsenBaseKernel/',dataName,'_CBKernel','.mat'],'KP');

labels = Y;
num_view = size(data,1);
num_samp = size(labels,1);
num_clus = length(unique(labels));


% initialize K_Lq
r = 12;

Hq = ones(num_view,1)/(2*num_view); 
% K_Lq = mycombFun(KP,Hq);  %% K_Lq is a matrix


variable1=[10];%参数范围
variable2=[10];%参数范围


knn = 20;    %%% obtain good results
normData = 1;
NITER = 30;


Max_Acc = 0;
% Max_NMI = 0;
% ACC_index = 0;

Numvariable1 = length(variable1);
Numvariable2 = length(variable2);

Numvariables = Numvariable1*Numvariable2;

% result = cell(Numvariable1,Numvariable2);

Result = cell(Numvariable1,Numvariable2);
OBJm = zeros(Numvariables,NITER);

for p=1:Numvariable1
    alpha = variable1(p);
    for q=1:Numvariable2
        beta = variable2(q);
        [result, Z, Tim, Wv, OBJ, H] = MVCsubspace_TM(data,labels, alpha, beta, knn, KP, Hq, NITER);
        
       disp(OBJ)
        
%         if size(OBJ,2) == size(OBJm,2)
% %             pq_OBJ = q+(p-1)*Numvariable1;
%             pq_OBJ = Numvariable1*Numvariable2;
%             OBJm(pq_OBJ,:) = OBJ;
%         end
        
%         OBJ(:,NITER) = OBJ(1,NITER);
        
%         for i= 1:Numvariables
%             OBJ(i,:) = OBJ;
%         end

        Acc = result(1);
        
        if Max_Acc < Acc
            Max_Acc = Acc;
%             ACC_index = ACC_index + 1;
%             Max_NMI = NMI;
            NMI = result(2);
            Purity = result(3);
            BesRes = [Max_Acc,NMI,Purity]; 
%         else
%             BesRes =  result;            
        end  
        Result{p,q} = result;
        pq_OBJ = q+(p-1)*Numvariable1;
        OBJm(pq_OBJ,:) = OBJ;
        TIM(pq_OBJ) = Tim;
        
        save([path,'myFinalRes/',dataName,'_Res-',num2str(alpha),'-',num2str(beta),'.mat'],'result','Tim');
        save([path,'myFinalRes/',dataName,'_OBJ','.mat'],'OBJm');     
    end
%     if Result{p,q} == BesRes
%         save([path,'BestRes/',dataName,'_Res-',num2str(alpha),'-',num2str(beta),'.mat'],'BesRes','Tim');
%     end
    
    
end
mean_TIM = mean(TIM);


% 
for p=1:Numvariable1
    alpha = variable1(p);
    for q=1:Numvariable2
        beta = variable2(q);
        if Result{p,q} == BesRes
           save([path,'BestRes/',dataName,'_Res-',num2str(alpha),'-',num2str(beta),'.mat'],'BesRes','mean_TIM');
           save([path,'BestRes/',dataName,'_Z'],'Z','H');
        end      
    end
end


