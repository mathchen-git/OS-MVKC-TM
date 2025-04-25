function Hq=solveHq(KP,num_view,Hq,Z)
% KPvZ=zeros(1,num_samp);
sumKPvZ=0;
for v=1:num_view        %计算核矩阵的和K
%     sumWK=sumWK+Ws(v)*Ks{v};%计算Ws*Ks的和
    KPvZ(v,1)=norm(Hq(v)*KP(v)-Z,'fro')^2;%计算Ks-K的F范数的平方
    sumKPvZ = sumKPvZ + KPvZ(v,1);
    
%     KPvZ(1,v)=norm(KP(v)-Z,'fro')^2;
%     KPvZ(v)=norm(KP(v)-Z,'fro')^2;
end


% Ws=(ones(1,r)-KsK/sumKsK)/(r-1);    %自己的Ws计算公式

 
% KPvZ(v)=norm(KP(v)-Z,'fro')^2;
 
% Hq(v)=(1-KPvZ(1,v)/sumKPvZ)/(num_view-1); 

Hq=(ones(num_view,1)-KPvZ/sumKPvZ)/(num_view-1); 
end



