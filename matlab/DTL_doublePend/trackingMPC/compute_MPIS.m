function compute_MPIS(Xc,Uc,Ak,K)
    % Xc,Uc: Polyherdons of constraints
    [F, G, ~] = convert_Poly2Mat(Xc, Uc);
    Fpi = @(i) (F+G*K)*Ak^i;
    Xpi = @(i) Polyhedron(Fpi(i), ones(size(Fpi(i), 1), 1));
    Xmpi = Xpi(0);
    i= 0;
    while(1) % 
        i = i + 1;
        Xmpi_tmp = and(Xmpi, Xpi(i));
        if Xmpi_tmp == Xmpi
            break;
        else
            Xmpi = Xmpi_tmp;
        end
    end
end
