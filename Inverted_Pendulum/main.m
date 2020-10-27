%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Model Predictive Control on Inverted Pendulum  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Discretization
dt = 0.01;	

% Inverted Pendulum Parameters
g = 9.81;       % gravity 
m = 2;          % Mass of the pendulum 
l = 1;          % length of pendulum
b = 0.1;        % damping coefficient
I = m*(l.^2);   % inertia 

% Initial Configuration: (Initial state)
xo = zeros(2,1);
xo(1,1) = 0;
xo(2,1) = 0;
current_state = xo; %Set current state to initial state
next_state = zeros(2,1);

% Target: (Terminal States)
p_target(1,1) = pi;     % theta
p_target(2,1) = 0;      % theta_dot

while current_state(1,1) <= 3.14
    u = fnDDP(current_state);
    %fprintf('Control Input %d \n',u(1,1));
    next_state = fnDynamics(current_state,u(1,1),dt); %Apply (only) first control input
    current_state = next_state;
    fprintf('Current State %d \n',current_state(1,1));
end
