function Z = updateZ(Z,num_samp,I,num_view,Wv,F,alpha,Zv)
    iniS = Z;
    Z = zeros(num_samp);
    A = 0*alpha*I;
    for v = 1:num_view
%         A = A + Wv(v)*L{v};   %%%  (15) and (16)
        A = A + alpha*Wv(v)*I;
    end
    dist_u = L2_distance_1(F',F');
    
    sumZ = zeros(num_samp);
    for v = 1:num_view
        sumZ = sumZ + Zv{v};
    end 
    
    for ni = 1:num_samp
        index = find(iniS(ni,:)>0);
%         b = 2*beta*I(ni,index) - lambda*dist_u(ni,index);
        b = 2*alpha*sumZ(ni,index) - dist_u(ni,index);
        % solve z^T*A*z-z^T*b
        [Zi, ~] = fun_alm(A(index,index),b);  %% update S
        Z(ni,index) = Zi';
    end
    Z = (Z+Z')/2; 
end

