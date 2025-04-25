%%%%%%% try to generate an index set 

clear
clc

X = [1 2 3 2; 3 2 1 1;3 4 1 0];
% B = 2*ones(3,4);
% X = A+B;
[m,n] = size(X);

sum_dXiXj = 0;
for i=1:n
    for j=1:n
        if j>=i
            dXiXj(i,j) = norm(X(:,i)-X(:,j),2).^2;
            sum_dXiXj = sum_dXiXj+dXiXj(i,j);
        end
    end   
end
sum_dXiXj2 = 2/(n*(n-1))*sum_dXiXj;

S = zeros(n,n);
for i=1:n
    for j=1:n   
        if dXiXj(i,j) > sum_dXiXj2
            S(i,j) = 1;
        else
            S(i,j) = 0;
        end
    end
end



