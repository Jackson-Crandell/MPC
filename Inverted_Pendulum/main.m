%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  MPC-DDP on Inverted Pendulum                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE8803  Fall  2020                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Jackson Crandell and Luis Pimental     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% Variables for the inverted pendulum
global g;
global m; 
global l; 
global I; 
global b;	

% Inverted Pendulum Parameters
g = 9.81;       % gravity 
m = 2;          % Mass of the pendulum 
l = 1;          % length of pendulum
b = 0.1;        % damping coefficient
I = m*(l.^2);   % inertia 

% Other Parameters
Horizon = 50;
num_iter = 20;
dt = 0.01;

% Weight in Final State: (part of terminal cost)
Q_f = zeros(2,2);
Q_f(1,1) = 100; 	%Penalize more w.r.t errors in theta
Q_f(2,2) = 10;      %Penalize more w.r.t errors in theta_dot

% Weight in the Control: 
R = .1* eye(1,1);

% Learning Rate .5
gamma = 0.5;

% Initial Configuration: (Initial state)
x_init = zeros(2,1);
x_init(1,1) = 0; %Initial Theta
x_init(2,1) = 0; %Initial Theta_dot
u_init = zeros(1,Horizon-1);

% Target: (Terminal States)
p_target(1,1) = pi;     % theta
p_target(2,1) = 0;      % theta_dot

%Final Trajectories
x_traj = zeros(2,Horizon);
u_traj = zeros(1,Horizon-1);
Cost = zeros(1,Horizon);

i = 0;
while 1
    [u, cost] = fnDDP(x_init,u_init,Horizon,num_iter,dt,p_target,gamma,Q_f,R);
    u_init = u; % Update controls
    [x] = fnsimulate(x_init,u,Horizon,dt,0); %Apply (only) first control input
    x_init = x(:,2);
    
    if (norm(x_init - p_target) < 1e-3) % Stop when target is reached
        break
    end
    
    if mod(i, 10) == 0
        fprintf('MPC Iteration %i,  Current Cost = %e \n',i,norm(x_init - p_target));
    end
    
    i = i + 1; % Keep track of number of iterations
    
    %Add to trajectories and cost
    x_traj(:,i) = x(:,2);
    u_traj(i) = u_init(1,1);
    Cost(i) = cost;
end

Horizon = length(x_traj);

time(1)=0;
for i= 2:Horizon
	time(i) =time(i-1) + dt;  
end

figure(1);
subplot(2,2,1);
hold on;
plot(time,x_traj(1,:),'linewidth',4);  
plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4);
title('$\theta$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20);
ylabel('Rad','fontsize',20);

hold off;
grid;

subplot(2,2,2);
hold on;
plot(time,x_traj(2,:),'linewidth',4); 
plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4);
title('$\dot{\theta}$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20);
ylabel('Rad/s','fontsize',20);

hold off;
grid;

subplot(2,2,3);
hold on
plot(Cost,'linewidth',2); 
xlabel('Iterations','fontsize',20);
title('Cost','fontsize',20);

subplot(2,2,4);
hold on;
title('Animation','fontsize',20);
invPend_animation(x_traj,Horizon);
hold off;
grid;
