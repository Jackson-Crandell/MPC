function [x] = fnsimulate(dynamicsf,xo,u_new,Horizon,dt,sigma)

x = xo;

for k = 1:(Horizon-1)
    x(:,k+1) = x(:,k) + dynamicsf(u_new(1,k),x(2,k),x(3,k),x(4,k))* dt; %* sqrt(dt) * sigma * randn ;
end