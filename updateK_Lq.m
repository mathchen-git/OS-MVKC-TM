% function K_Lq = updateK_Lq(num_samp,I,num_view,Z,Zv,alpha)
% function K_Lq = updateK_Lq(num_samp,num_view,Z,Zv)
function KP2 = updateK_Lq(num_samp,num_view,Z,Zv,alpha,Hq,I)
    sumZv = zeros(num_samp);
    sumZvZvT = zeros(num_samp);

for v = 1:num_view
    KP2{v} = 0.5*(1/(alpha*Hq(v)*Hq(v))+0.0001)*(2*alpha*Hq(v)*Z+2*Zv{v}-2*Zv{v}*Zv{v}'-I);
    KP2{v} = (KP2{v}+KP2{v}')/2;
end    
    
end

