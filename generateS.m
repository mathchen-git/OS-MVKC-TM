function S = generateS(Z)
n = size(Z,1);

sum_dZiZj = 0;
for i=1:n
    for j=1:n
        if j>=i
            dZiZj(i,j) = norm(Z(:,i)-Z(:,j),2).^2;
            sum_dZiZj = sum_dZiZj+dZiZj(i,j);
        end
    end   
end
sum_dZiZj2 = 2/(n*(n-1))*sum_dZiZj;

S = zeros(n,n);
for i=1:n
    for j=1:n   
        if dZiZj(i,j) > sum_dZiZj2
            S(i,j) = 1;
        else
            S(i,j) = 0;
        end
    end
end
end

