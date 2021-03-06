%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  MPC-DDP on Inverted Pendulum                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE8803  Fall  2020                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Jackson Crandell and Luis Pimental     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% Variables for the inverted pendulum
global g_nominal;
global m_nominal; 
global l_nominal; 
global I_nominal; 
global b_nominal;	

% Inverted Pendulum Parameters
g_nominal = 9.81;       % gravity 
m_nominal = 2;          % Mass of the pendulum 
l_nominal = 1;          % length of pendulum
b_nominal = 0.1;        % damping coefficient
I_nominal = m_nominal*(l_nominal.^2);   % inertia 

global g_real; 
global m_real;
global l_real; 
global b_real; 
global I_real; 

sigma = 0.0; 

sigma_levels = [];

real_g_params = [];
real_m_params = [];
real_b_params = [];
real_l_params = [];
real_I_params = [];   



while 1
    
    g_real = abs(g_nominal + g_nominal*(0.5^(sigma))*(rand*2 - 1)*sigma^2);  		
    m_real = abs(m_nominal + m_nominal*(rand*2 - 1)*sigma);            
    l_real = abs(l_nominal + l_nominal*(rand*2 - 1)*sigma);  
    b_real = abs(b_nominal + b_nominal*(rand*2 - 1)*sigma);        
    I_real = abs(I_nominal + I_nominal*(rand*2 - 1)*sigma);     
    
    syms x1 x2 u

    F_nominal = [x2, 
        ((-b_nominal/I_nominal)*x2)-((m_nominal*g_nominal*l_nominal)/I_nominal)*sin(x1)+(u/I_nominal)
    ]; 
    dynamicsf_nominal = matlabFunction(F_nominal);


    F_real = [x2, 
        ((-b_real/I_real)*x2)-((m_real*g_real*l_real)/I_real)*sin(x1)+(u/I_real)
    ];
    dynamicsf_real = matlabFunction(F_real);
    
    
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

    iter = 0;
    while 1
        [u, cost] = fnDDP(dynamicsf_nominal,x_init,u_init,Horizon,num_iter,dt,p_target,gamma,Q_f,R);
        u_init = u; % Update controls
        [x] = fnsimulate(dynamicsf_real,x_init,u,Horizon,dt,0); %Apply (only) first control input
        x_init = x(:,2);

        if (norm(x_init - p_target) <= 1e-3) % Stop when target is reached
            break
        end
        
        if(iter == 500 && norm(x_init - p_target) > 1e-1)
            break   
        end
        
        
        if mod(iter, 10) == 0
            fprintf('MPC Iteration %i, Current Cost = %e \n',iter,norm(x_init - p_target));
        end

        iter = iter + 1; % Keep track of number of iterations

        %Add to trajectories and cost
        x_traj(:,iter) = x(:,2);
        u_traj(iter) = u_init(1,1);
        Cost(iter) = cost;
    end
    
    Horizon = length(x_traj);

    time = NaN
    time(1:Horizon) = 0;

    for i= 2:Horizon
        time(i) =time(i-1) + dt;  
    end
    
    
    figure(1);
    subplot(2,2,1);
    hold on;
    plot(time,x_traj(1,:),'linewidth',4);  
    title('$\theta$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20);
    ylabel('Rad','fontsize',20);
    grid;

    subplot(2,2,2);
    hold on;
    plot(time,x_traj(2,:),'linewidth',4); 
    title('$\dot{\theta}$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20);
    ylabel('Rad/s','fontsize',20);
    grid;

    subplot(2,2,3);
    hold on
    plot(Cost,'linewidth',2); 
    xlabel('Iterations','fontsize',20);
    title('Cost','fontsize',20);
    
    sigma_levels = [sigma_levels sigma]
    real_g_params = [real_g_params g_real]
    real_m_params = [real_m_params m_real]
    real_b_params = [real_b_params b_real]
    real_l_params = [real_l_params l_real]
    real_I_params = [real_I_params I_real]

    sigma = sigma + 0.4;
    iter 
    
    if(iter == 500)
        break
    end
    
end 

%%

figure(1);
subplot(2,2,1)
hold on
plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)

subplot(2,2,2);
hold on
plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)


figure(1) 
subplot(2,2,1)
legend('$\sigma = 0$','$\sigma = 0.4$', '$\sigma = 0.8$','$\sigma = 1.2$', '$\sigma = 1.6$','$\sigma = 2.0$','$\sigma = 2.4$','$\sigma = 2.8$','Target',  'Interpreter','latex','fontsize',24)

subplot(2,2,2)
legend('$\sigma = 0$','$\sigma = 0.4$', '$\sigma = 0.8$','$\sigma = 1.2$', '$\sigma = 1.6$','$\sigma = 2.0$','$\sigma = 2.4$','$\sigma = 2.8$','Target',  'Interpreter','latex','fontsize',24)

%%


figure(2);
subplot(3,3,1)
hold on
plot(time,g_nominal*ones(1,Horizon),'red','linewidth',4); 
title('$g$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('$m/s^2$','Interpreter','latex','fontsize',20);
grid;

subplot(3,3,2)
hold on
plot(time,m_nominal*ones(1,Horizon),'red','linewidth',4); 
title('$m$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('$Kg$','Interpreter','latex','fontsize',20);
grid;

subplot(3,3,3)
hold on
plot(time,b_nominal*ones(1,Horizon),'red','linewidth',4); 
title('$b$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('$Kg$','Interpreter','latex','fontsize',20);
hold off;
grid;

subplot(3,3,4)
hold on
plot(time,l_nominal*ones(1,Horizon),'red','linewidth',4); 
title('$l$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('Meters','fontsize',20);
hold off;
grid;

subplot(3,3,5)
hold on
plot(time,I_nominal*ones(1,Horizon),'red','linewidth',4); 
title('$l$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('Meters','fontsize',20);
hold off;
grid;

for i = 1:length(sigma_levels)
    g_plot = real_g_params(i);
    m_plot = real_m_params(i);
    b_plot = real_b_params(i);
    l_plot = real_l_params(i);
    I_plot = real_I_params(i);
    
    figure(2);
    subplot(3,3,1)
    hold on
    plot(time,g_plot*ones(1,Horizon),'linewidth',4); 
    title('$g$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('$m/s^2$','Interpreter','latex','fontsize',20);
    grid;
    
    subplot(3,3,2)
    hold on
    plot(time,m_plot*ones(1,Horizon),'linewidth',4); 
    title('$m$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('$Kg$','Interpreter','latex','fontsize',20);
    grid;
    
    subplot(3,3,3)
    hold on
    plot(time,b_plot*ones(1,Horizon),'linewidth',4); 
    title('$b$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('$Kg$','Interpreter','latex','fontsize',20);
    hold off;
    grid;
    
    subplot(3,3,4)
    hold on
    plot(time,l_plot*ones(1,Horizon),'linewidth',4); 
    title('$l$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('Meters','fontsize',20);
    hold off;
    grid;
    
    subplot(3,3,5)
    hold on
    plot(time,I_plot*ones(1,Horizon),'linewidth',4); 
    title('$I$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('Meters','fontsize',20);
    hold off;
    grid;
    
end


figure(2) 
subplot(3,3,1)
legend('Nominal','$\sigma = 0$','$\sigma = 0.4$', '$\sigma = 0.8$','$\sigma = 1.2$', '$\sigma = 1.6$','$\sigma = 2.0$','$\sigma = 2.4$','$\sigma = 2.8$','Interpreter','latex','fontsize',24)

subplot(3,3,2)
legend('Nominal','$\sigma = 0$','$\sigma = 0.4$', '$\sigma = 0.8$','$\sigma = 1.2$', '$\sigma = 1.6$','$\sigma = 2.0$','$\sigma = 2.4$','$\sigma = 2.8$','Interpreter','latex','fontsize',24)

subplot(3,3,3)
legend('Nominal','$\sigma = 0$','$\sigma = 0.4$', '$\sigma = 0.8$','$\sigma = 1.2$', '$\sigma = 1.6$','$\sigma = 2.0$','$\sigma = 2.4$','$\sigma = 2.8$','Interpreter','latex','fontsize',24)

subplot(3,3,4)
legend('Nominal','$\sigma = 0$','$\sigma = 0.4$', '$\sigma = 0.8$','$\sigma = 1.2$', '$\sigma = 1.6$','$\sigma = 2.0$','$\sigma = 2.4$','$\sigma = 2.8$','Interpreter','latex','fontsize',24)

subplot(3,3,5)
legend('Nominal','$\sigma = 0$','$\sigma = 0.4$', '$\sigma = 0.8$','$\sigma = 1.2$', '$\sigma = 1.6$','$\sigma = 2.0$','$\sigma = 2.4$','$\sigma = 2.8$','Interpreter','latex','fontsize',24)



















