function [Xmpi]=compute_initS(Xc,Ak)
    % Xc: Polyherdon of constraints
    % Ak: Matrix of the autonomous system under Xc constraints
    
    % normalize lhs w.r.t. rhs vector in H-representation
    rhe2one = @(poly) poly.A./repmat(poly.b, 1, size(poly.A, 2));
    F=rhe2one(Xc);
    
    Fpi = @(i) (F)*Ak^i;
    Xpi = @(i) Polyhedron(Fpi(i), ones(size(Fpi(i), 1), 1));
    Xmpi = Xpi(0);
    i= 0;
    while(i<2) % 
        i = i + 1;
        % compute set intersection
        Xmpi_tmp = Xmpi & Xpi(i);
        if Xmpi_tmp == Xmpi
            break;
        else
            Xmpi = Xmpi_tmp;
        end
    end
    fprintf("Interations : %d\n",i);
end
