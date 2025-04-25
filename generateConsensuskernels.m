

function KP=generateConsensuskernels(KH1G,KH2G,KH3G,KH4G,KH5G,KH6G,KH7G,KHC,KH1M,KH2M,KH3M,KH4M,Lq)  %%%%%generate consensus base kernels K1, K2, ..., Kp, ..., Kr-------generatemultiplekernel
[n,d,numView]=size(KH1G); %%% KH1G1 is a matrix with n*n*view
% Kbeta=zeros(n,d);
KP=zeros(n,d,numView);

for p=1:numView

%         KM(:,:,i)=KM(:,:,i)+alpha(i,j)*KHG_KHM_KHC{i}(:,:,j);
%         KM(:,:,p)=alpha(p,1)*KHG(:,:,p)+alpha(p,2)*KHM(:,:,p)+alpha(p,3)*KHC(:,:,p); %% KM:n*n*view  
        KP(:,:,p)=Lq(p,1)*KH1G(:,:,p)+Lq(p,2)*KH2G(:,:,p)+Lq(p,3)*KH3G(:,:,p)+Lq(p,4)*KH4G(:,:,p)...
                  +Lq(p,5)*KH5G(:,:,p)+Lq(p,6)*KH6G(:,:,p)++Lq(p,7)*KH7G(:,:,p)+Lq(p,8)*KHC(:,:,p)...
                  +Lq(p,9)*KH1M(:,:,p)+Lq(p,10)*KH2M(:,:,p)+Lq(p,11)*KH3M(:,:,p)+Lq(p,12)*KH4M(:,:,p); %% KM:n*n*view 
end



%% �Ҿ�����������Ļ��ͱ�׼����������Ҫ�ģ����ǲ�ȷ��
%KM = kcenter(KM);
%KM = knorm(KM);

%% ������飬myconbFun���ˣ���mkkmeans_train��ʹ����

%Kbeta= kcenter(Kbeta);
%Kbeta = knorm(Kbeta);