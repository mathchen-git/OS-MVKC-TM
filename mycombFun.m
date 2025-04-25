function K_Lq = mycombFun(KP,Hq)   %%% construct K_beta

m = size(KP,3);
n = size(KP,1);
K_Lq = zeros(n);
for p =1:m
    K_Lq = K_Lq + KP(:,:,p)*Hq(p);   %%% gamma is beta in paper
end