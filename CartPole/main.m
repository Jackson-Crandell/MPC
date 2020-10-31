%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  MPC-DDP on CartPole                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE4803  Fall  2020                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Jackson Crandell and Luis Pimental     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

save_video = false;

% Variables for the cart-pole problem
global g; 
global m_c;
global m_p; 
global l; 

g = 9.81;  		% gravity 
m_c = 1.0; 		% mass of the cart
m_p = 0.01;		% mass of the pole
l = 0.25;  		% length of the pole

%Dynamics
% Initialize the state space equations for the dynamics 
% x1: dot_x (velocity)
% x2: ddot_x (acceleration)
% x3: dot_theta (angular velocity)
% x4: ddot_theta (angular accelaration)
% u: control input (force)
syms x1 x2 x3 x4 u

F = [  x2,
        (m_p*sin(x3)*(l*(x4.^2) + g*cos(x3)) + u)./(m_c + m_p*(sin(x3).^2)), 
        x4, 
        (-m_p*l*(x4.^2)*cos(x3)*sin(x3) - (m_c + m_p)*g*sin(x3) - cos(x3)*u)./(l*(m_c + m_p*(sin(x3).^2)))
]; 

dynamicsf = matlabFunction(F);


% Horizon 
Horizon = 60; 

% Number of Iterations
num_iter = 30;

% Discretization
dt = 0.01;

% Weight in Final State: (part of terminal cost)
Q_f = eye(4,4);
Q_f(1,1) = 100;		% Penalize more w.r.t errors in position
Q_f(2,2) = 100;     % Penalize less for errors in velocity
Q_f(3,3) = 10;    	% Penalize less for errors in anglular position
Q_f(4,4) = 10;      % Penalize less for errors in angular velocity

% Weight in the Control:
R = 1* eye(1,1); 	

% Initial Configuration: (Initial state)
x_init = zeros(4,1);
x_init(1,1) = 0;
x_init(2,1) = 0;
x_init(3,1) = 0;
x_init(4,1) = 0;
u_init = zeros(1,Horizon-1);

% Target: (Terminal States)
p_target(1,1) = 0;          % x         (position of cart)
p_target(2,1) = 0;          % x_dot     (velocity of the cart)
p_target(3,1) = pi;         % theta     (angular position of pole)
p_target(4,1) = 0;          % theta_dot (angular velocity of the pole)

% Learning Rate
gamma = 0.4; 

%Final Trajectories
x_traj = zeros(4,Horizon);
u_traj = zeros(2,Horizon-1);
Cost = zeros(1,Horizon);
i = 0;

while 1
    [u_new, cost] = fnDDP(dynamicsf,x_init,u_init,num_iter,Horizon,gamma,p_target,dt,Q_f,R);
    u_init = u_new;
    [x] = fnsimulate(dynamicsf,x_init,u_new,Horizon,dt,0);
    x_init = x(:,2);
    
    if(norm(x_init - p_target) < 1e-3) %Stop when task is complete
        break
    end
    i = i + 1; % Keep track of iterations
    
    if mod(i, 10) == 0
        fprintf('MPC Iteration %i,  Current Cost = %e \n',i,norm(x_init - p_target));
    end 

    %Add to trajectories and cost
    x_traj(:,i) = x(:,4);
    u_traj(1,i) = u_init(1,1);
    Cost(1,i) = cost;
end

Horizon = length(x_traj);

time(1)=0;
for i= 2:Horizon
	time(i) =time(i-1) + dt;  
end

for k=1:5:Horizon
    drawCartpend(x_traj(:,k),m_p,m_c,l);
    if (save_video)
        E(k) = getframe(gcf);
    end
end

% Save Video
if (save_video)
    E(1).cdata = E(1).cdata(:,1:1000,:); %Resize Array to make video
    video = VideoWriter('CartPole.avi','Uncompressed AVI');
    open(video)
    writeVideo(video,E)
    close(video)
end

figure(1);
subplot(3,2,1)
hold on
plot(time,x_traj(1,:),'linewidth',4);  
plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)
title('$x$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('Meters','fontsize',20);
hold off;
grid;


subplot(3,2,2);hold on;
plot(time,x_traj(2,:),'linewidth',4); 
plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)
title('$\dot{x}$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('Meters/s','fontsize',20);
hold off;
grid;

subplot(3,2,3);hold on
plot(time,x_traj(3,:),'linewidth',4); 
plot(time,p_target(3,1)*ones(1,Horizon),'red','linewidth',4)
title('$\theta$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('Rad','fontsize',20);
hold off;
grid;

subplot(3,2,4);hold on
plot(time,x_traj(4,:),'linewidth',4); 
plot(time,p_target(4,1)*ones(1,Horizon),'red','linewidth',4)
title('$\dot{\theta}$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('Rad/s','fontsizecch',20);
hold off;
grid;

subplot(3,2,5);hold on
plot(Cost,'linewidth',2); 
xlabel('Iterations','fontsize',20)
title('Cost','fontsize',20);
