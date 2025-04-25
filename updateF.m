
function [F,ev,F_old] = updateF(S,F,num_clus)

LS = diag(sum(S)) - S;
F_old = F;
[F, ~, ev] = eig1(LS, num_clus, 0);

end
