function [A,B] = Jacobians(x,u, dFx, dFu)
    A = dFx(u(1,1),x(3,1), x(4,1)); 
    B = dFu(x(3,1));
end
