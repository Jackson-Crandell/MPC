%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  MPC-DDP on CartPole                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE4803  Fall  2020                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Jackson Crandell and Luis Pimental     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_video = false;

% Variables for the cart-pole problem
global g_nominal; 
global m_c_nominal;
global m_p_nominal; 
global l_nominal; 

g_nominal = 9.81;  		% gravity 
m_c_nominal = 1.0; 		% mass of the cart
m_p_nominal = 0.01;		% mass of the pole
l_nominal = 0.25;  		% length of the pole

global g_real; 
global m_c_real;
global m_p_real; 
global l_real; 

sigma = 0.0; 
sigma_levels = []
real_g_params = []
real_m_c_params = []
real_m_p_params = []
real_l_params = []

while(1)
    g_real = g_nominal + g_nominal*0.5*(rand*2 - 1)*sigma^2;  		
    m_c_real = m_c_nominal + m_c_nominal*(rand*2 - 1)*sigma;            
    m_p_real = m_p_nominal + m_p_nominal*(rand*2 - 1)*sigma;  
    l_real = l_nominal + l_nominal*(rand*2 - 1)*sigma;        
    


    %Dynamics
    % Initialize the state space equations for the dynamics 
    % x1: dot_x (velocity)
    % x2: ddot_x (acceleration)
    % x3: dot_theta (angular velocity)
    % x4: ddot_theta (angular accelaration)
    % u: control input (force)

    syms x1 x2 x3 x4 u

    F_nominal = [  x2,
            (m_p_nominal*sin(x3)*(l_nominal*(x4.^2) + g_nominal*cos(x3)) + u)./(m_c_nominal + m_p_nominal*(sin(x3).^2)), 
            x4, 
            (-m_p_nominal*l_nominal*(x4.^2)*cos(x3)*sin(x3) - (m_c_nominal + m_p_nominal)*g_nominal*sin(x3) - cos(x3)*u)./(l_nominal*(m_c_nominal + m_p_nominal*(sin(x3).^2)))
    ]; 

    dynamicsf_nominal = matlabFunction(F_nominal);

    F_real = [  x2,
            (m_p_real*sin(x3)*(l_real*(x4.^2) + g_real*cos(x3)) + u)./(m_c_real + m_p_real*(sin(x3).^2)), 
            x4, 
            (-m_p_real*l_real*(x4.^2)*cos(x3)*sin(x3) - (m_c_real + m_p_real)*g_real*sin(x3) - cos(x3)*u)./(l_real*(m_c_real + m_p_real*(sin(x3).^2)))
    ]; 

    dynamicsf_real = matlabFunction(F_real);


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
    iter = 0;
   
    while 1
        [u_new, cost] = fnDDP(dynamicsf_nominal,x_init,u_init,num_iter,Horizon,gamma,p_target,dt,Q_f,R);
        u_init = u_new;
        [x] = fnsimulate(dynamicsf_real, x_init,u_new,Horizon,dt,0);
        x_init = x(:,2);

        if(norm(x_init - p_target) <= 1e-2) %Stop when task is complete
            break
        end

        if(iter == 1000 && norm(x_init - p_target) > 1e-2)
            break   
        end
        
        iter = iter + 1; % Keep track of iterations

        if mod(iter, 10) == 0
            fprintf('MPC Iteration %i,  Current Cost = %e \n',iter,norm(x_init - p_target));
        end 

        %Add to trajectories and cost
        x_traj(:,iter) = x(:,2);
        u_traj(1,iter) = u_init(1,1);
        Cost(1,iter) = cost;
    end

    Horizon = length(x_traj);

    time = NaN
    time(1:Horizon) = 0;

    for i= 2:Horizon
        time(i) =time(i-1) + dt;  
    end

    figure(1);
    subplot(3,2,1)
    hold on
    plot(time,x_traj(1,:),'linewidth',4);  
    title('$x$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('Meters','fontsize',20);
    hold off;
    grid;


    subplot(3,2,2);hold on;
    plot(time,x_traj(2,:),'linewidth',4); 
    hold on
    title('$\dot{x}$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('Meters/s','fontsize',20);
    hold off;
    grid;

    subplot(3,2,3);hold on
    plot(time,x_traj(3,:),'linewidth',4); 
    hold on
    title('$\theta$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('Rad','fontsize',20);
    hold off;
    grid;

    subplot(3,2,4);hold on
    plot(time,x_traj(4,:),'linewidth',4); 
    hold on
    title('$\dot{\theta}$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('Rad/s','fontsize',20);
    hold off;
    grid;

    subplot(3,2,5);
    hold on
    plot(Cost,'linewidth',2); 
    xlabel('Iterations','fontsize',20)
    title('Cost','fontsize',20);

    sigma_levels = [sigma_levels sigma]
    real_g_params = [real_g_params g_real]
    real_m_c_params = [real_m_c_params m_c_real]
    real_m_p_params = [real_m_p_params m_p_real]
    real_l_params = [real_l_params l_real]

    sigma = sigma + 0.2;

    if(iter == 1000)
        break
    end

end



%%

figure(1);
subplot(3,2,1)
hold on
plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)

subplot(3,2,2);
hold on
plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)

subplot(3,2,3);
hold on
plot(time,p_target(3,1)*ones(1,Horizon),'red','linewidth',4)

subplot(3,2,4);
hold on
plot(time,p_target(4,1)*ones(1,Horizon),'red','linewidth',4)

subplot(3,2,5);
hold on
plot(Cost,'linewidth',2); 
xlabel('Iterations','fontsize',20)
title('Cost','fontsize',20);

figure(1) 
subplot(3,2,1)
legend('$\sigma = 0$','$\sigma = 0.2$', '$\sigma = 0.4$','$\sigma = 0.6$', '$\sigma = 0.8$','$\sigma = 1.0$','$\sigma = 1.2$','Target',  'Interpreter','latex','fontsize',24)

subplot(3,2,2)
legend('$\sigma = 0$','$\sigma = 0.2$', '$\sigma = 0.4$','$\sigma = 0.6$', '$\sigma = 0.8$','$\sigma = 1.0$','$\sigma = 1.2$','Target',  'Interpreter','latex','fontsize',24)

subplot(3,2,3)
legend('$\sigma = 0$','$\sigma = 0.2$', '$\sigma = 0.4$','$\sigma = 0.6$', '$\sigma = 0.8$','$\sigma = 1.0$','$\sigma = 1.2$','Target',  'Interpreter','latex','fontsize',24)

subplot(3,2,4)
legend('$\sigma = 0$','$\sigma = 0.2$', '$\sigma = 0.4$','$\sigma = 0.6$', '$\sigma = 0.8$','$\sigma = 1.0$','$\sigma = 1.2$','Target',  'Interpreter','latex','fontsize',24)

subplot(3,2,5)
legend('$\sigma = 0$','$\sigma = 0.2$', '$\sigma = 0.4$','$\sigma = 0.6$', '$\sigma = 0.8$','$\sigma = 1.0$','$\sigma = 1.2$',  'Interpreter','latex','fontsize',24)

%%
figure(2);
subplot(2,2,1)
hold on
plot(time,g_nominal*ones(1,Horizon),'red','linewidth',4); 
title('$g$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('$m/s^2$','Interpreter','latex','fontsize',20);
grid;

subplot(2,2,2)
hold on
plot(time,m_c_nominal*ones(1,Horizon),'red','linewidth',4); 
title('$m_p$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('$Kg$','Interpreter','latex','fontsize',20);
grid;

subplot(2,2,3)
hold on
plot(time,m_p_nominal*ones(1,Horizon),'red','linewidth',4); 
title('$m_p$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('$Kg$','Interpreter','latex','fontsize',20);
hold off;
grid;

subplot(2,2,4)
hold on
plot(time,l_nominal*ones(1,Horizon),'red','linewidth',4); 
title('$l$','Interpreter','latex','fontsize',24);
xlabel('Time in sec','fontsize',20)
ylabel('Meters','fontsize',20);
hold off;
grid;

for i = 1:length(sigma_levels)
    g_plot = real_g_params(i);
    m_c_plot = real_m_p_params(i);
    m_p_plot = real_m_p_params(i);
    l_plot = real_l_params(i);
    
    figure(2);
    subplot(2,2,1)
    hold on
    plot(time,g_plot*ones(1,Horizon),'linewidth',4); 
    title('$g$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('$m/s^2$','Interpreter','latex','fontsize',20);
    grid;
    
    subplot(2,2,2)
    hold on
    plot(time,m_c_plot*ones(1,Horizon),'linewidth',4); 
    title('$m_p$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('$Kg$','Interpreter','latex','fontsize',20);
    grid;
    
    subplot(2,2,3)
    hold on
    plot(time,m_p_plot*ones(1,Horizon),'linewidth',4); 
    title('$m_p$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('$Kg$','Interpreter','latex','fontsize',20);
    hold off;
    grid;
    
    subplot(2,2,4)
    hold on
    plot(time,l_plot*ones(1,Horizon),'linewidth',4); 
    title('$l$','Interpreter','latex','fontsize',24);
    xlabel('Time in sec','fontsize',20)
    ylabel('Meters','fontsize',20);
    hold off;
    grid;
    
end


figure(2) 
subplot(2,2,1)
legend('Nominal','$\sigma = 0$','$\sigma = 0.2$', '$\sigma = 0.4$','$\sigma = 0.6$', '$\sigma = 0.8$','$\sigma = 1.0$','$\sigma = 1.2$','Interpreter','latex','fontsize',24)

subplot(2,2,2)
legend('Nominal','$\sigma = 0$','$\sigma = 0.2$', '$\sigma = 0.4$','$\sigma = 0.6$', '$\sigma = 0.8$','$\sigma = 1.0$','$\sigma = 1.2$','Interpreter','latex','fontsize',24)

subplot(2,2,3)
legend('Nominal','$\sigma = 0$','$\sigma = 0.2$', '$\sigma = 0.4$','$\sigma = 0.6$', '$\sigma = 0.8$','$\sigma = 1.0$','$\sigma = 1.2$','Interpreter','latex','fontsize',24)

subplot(2,2,4)
legend('Nominal','$\sigma = 0$','$\sigma = 0.2$', '$\sigma = 0.4$','$\sigma = 0.6$', '$\sigma = 0.8$','$\sigma = 1.0$','$\sigma = 1.2$','Interpreter','latex','fontsize',24)



