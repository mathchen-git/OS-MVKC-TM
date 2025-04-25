% function Wv = updateWv(num_view,K_Lq,L)
function Wv = updateWv(num_view,Z,L)
    for v = 1:num_view
        Wv(v) = 1/(2*sqrt(trace(Z'*L{v}*Z))+0.0001); % (0,1]
    end
end

