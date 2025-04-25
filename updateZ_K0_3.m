% function Z = updateZ_K0_3(Z,num_samp,I,num_view,Wv,L,F,alpha,beta,K_Lq)
function Z = updateZ_K0_3(Z,num_samp,I,num_view,Wv,L,F,alpha,beta,KP,Hq)
    iniZ = Z;
    Z = zeros(num_samp);
%     A = beta*I;
    A_2 = alpha*num_view*I;
    sumlvKi=zeros(num_samp);
    for v = 1:num_view
        A_1 = Wv(v)*L{v};   %%%  (15) and (16)
        sumlvKi=sumlvKi+Hq(v)*KP(v);
    end
    A = A_1+A_2; 
    dist_u = L2_distance_1(F',F');
    for ni = 1:num_samp
        index = find(iniZ(ni,:)>0);
        b =  2*alpha*sumlvKi(ni,index)- beta*dist_u(ni,index);
        % solve z^T*A*z-z^T*b
        [si, ~] = fun_alm(A(index,index),b);  %% update S
        Z(ni,index) = si';
    end
    Z = (Z+Z')/2; 
end

