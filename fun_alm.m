function [v, obj] = fun_alm(A,b)
if size(b,1) == 1
    b = b';
end

% initialize
rho = 1.5;
mu = 30;
n = size(A,1);
alpha = ones(n,1);
v = ones(n,1)/n;
% obj_old = v'*A*v-v'*b;

obj = [v'*A*v-v'*b];  % (16)
iter = 0;
while iter < 10
    % update z
    z = v-A'*v/mu+alpha/mu;  %  (20), z--p, v--s_i, mu--eta, alpha--q
    
    % update v
    c = A*z-b;        % A*p-b
    d = alpha/mu-z;   % (1/eta)*q-p
    mm = d+c/mu;      % (21)
    v = EProjSimplex_new(-mm);   % solve (21), v--s_i 
    
    % update alpha and mu
    alpha = alpha+mu*(v-z);   % alpha--q, mu--eta, v--s_i, z--p
    mu = rho*mu;    %  mu--eta,
    iter = iter+1;
    obj = [obj;v'*A*v-v'*b];
end
end

