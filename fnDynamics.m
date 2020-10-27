function next_state = fnDynamics(xo,u_new,dt)

global g;
global m; 
global l; 
global I; 
global b; 

x = xo;
Fx(1,1) = x(2,1); 
Fx(2,1) = ((-b/I)*x(2,1))-((m*g*l)/I)*sin(x(1,1));
G_x(2,1) = 1/I;

%Return only next state
next_state = x(:,1) + Fx * dt + G_x * u_new(1,:) * dt; % + G_x * u_new(:,1) * sqrt(dt) * sigma * randn; %Creates trajectory based on dynamics