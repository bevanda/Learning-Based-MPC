function [Xmpi]=compute_MPIS(Xw,Aw)
    % Xw: MPT Polyherdon object of constraints
    % Aw: Matrix of the autonomous system under Xw constraints
    
    % normalize lhs w.r.t. rhs vector in H-representation
    % to have  { x | F*x <= 1 }
    rhe2one = @(poly) poly.A./repmat(poly.b, 1, size(poly.A, 2));
    F=rhe2one(Xw);
    Fpi = @(i) (F)*Aw^i;
    Xpi = @(i) Polyhedron(Fpi(i), ones(size(Fpi(i), 1), 1));
    Xmpi = Xpi(0);
    i= 0;
    while(1) % 
        i = i + 1;
        % compute set intersection
        Xmpi_tmp = Xmpi & Xpi(i);
        if Xmpi_tmp == Xmpi
            break;
        else
            Xmpi = Xmpi_tmp;
        end
    end
    fprintf("Iterations : %d\n",i);
end
