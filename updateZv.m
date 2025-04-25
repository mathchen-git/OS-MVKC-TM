function [Zv,L_Zv] = updateZv(num_view,Wv,KP,Dv,Z,I)

    Zv = cell(num_view,1);
    L_Zv = cell(num_view,1);
  
    
    for v = 1:num_view
        G = (Dv{v})^(-0.5)*Z;
        Zv{v} = (2*KP(v)+0.0001*I)\(2*KP(v)+Wv(v)*(G*G'));
        Zv{v}(Zv{v}<0) = 0; 
        Zv{v} = (Zv{v} + Zv{v}')/2;
        Zv{v} = Zv{v} - diag(diag(Zv{v}));
        L_Zv{v} = diag(sum(Zv{v})) - Zv{v};
    end
end

