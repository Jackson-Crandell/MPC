function [A,B] = Jacobians(x,u,dFx,dFu)
    A = dFx(x(1,1)); 
    B = dFu();
end 
