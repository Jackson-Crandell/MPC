function  [u_traj] = fnDDP(current_state)

Horizon = 300; % 1.5sec

% Number of Iterations
num_iter = 100;

% Discretization
dt = 0.01;	

% Weight in Final State: (part of terminal cost)
Q_f = zeros(2,2);
Q_f(1,1) = 100; 	%Penalize more w.r.t errors in theta
Q_f(2,2) = 10;      %Penalize more w.r.t errors in theta_dot

%  Weight in the state (for running cost)
P = zeros(2,2); 
P(1,1) = 0;
P(2,2) = 0;

% Weight in the Control: 
R = .1* eye(1,1); 

% Initial Configuration: (Initial state)
xo = zeros(2,1);
xo(1,1) = current_state(1,1);
xo(2,1) = current_state(2,1);

% Initial Control:
% Horizon -1 b/c at last time step we have no control
u_k = zeros(1,Horizon-1);   
du_k = zeros(1,Horizon-1);

% Initial trajectory:
x_traj = zeros(2,Horizon);

% Target: (Terminal States)
p_target(1,1) = pi;     % theta
p_target(2,1) = 0;      % theta_dot

% Learning Rate .5
gamma = 0.5;

%Initialize Q Value Function
Q = zeros(1,Horizon);
Q_x = zeros(2,Horizon);
Q_u = zeros(1,Horizon);
Q_xx = zeros(2,2,Horizon);
Q_uu = zeros(1,1,Horizon);
Q_ux = zeros(1,2,Horizon);
 
for k = 1:num_iter % Run for a certain number of iterations

    %------------------------------------------------> Linearization of the dynamics
    %------------------------------------------------> Quadratic Approximations of the cost function 
    for  j = 1:(Horizon-1) %Discretize trajectory for each timestep

        % Linearization of the dynamics
        [dfx, dfu] = Jacobians(x_traj(:,j),u_k(:,j));

        A(:,:,j) = eye(2) + dfx * dt;    
        B(:,:,j) = dfu * dt;     

        % Quadratic expansion of the running cost around the x_traj (nominal trajectory) and u_k (nominal control)
        [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x_traj(:,j), u_k(:,j), j,R,dt); 

        L0(j) = dt * l0;            
        Lx(:,j) = dt * l_x;         
        Lxx(:,:,j) = dt * l_xx;     

        Lu(:,j) = dt * l_u;         
        Luu(:,:,j) = dt * l_uu;     
        Lux(:,:,j) = dt * l_ux; 

    end

    %------------------------------------------------> Boundary Conditions
    % Initialize value function at the boundary conditions
    Vxx(:,:,Horizon)= Q_f;                                                                                                                                                               
    Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target);                                       
    V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target); 


    %------------------------------------------------> Backpropagation of the Value Function
    for j = (Horizon-1):-1:1

         Q = L0(j) + V(:,j+1);
         Q_x = Lx(:,j) + A(:,:,j)'*Vx(:,j+1);
         Q_xx = Lxx(:,:,j) + A(:,:,j)'*Vxx(:,:,j+1)*A(:,:,j);
         Q_u  = Lu(:,j) + B(:,:,j)'*Vx(:,j+1);
         Q_uu = Luu(:,:,j) + B(:,:,j)'*Vxx(:,:,j+1)*B(:,:,j);
         Q_ux = Lux(:,:,j) + B(:,:,j)'*Vxx(:,:,j+1)*A(:,:,j);

         inv_Q_uu = inv(Q_uu);
         L_k(:,:,j)= - inv_Q_uu*Q_ux;   % Feedback term
         l_k (:,j) = - inv_Q_uu*Q_u;    % Feedforward term

         Vxx(:,:,j) = Q_xx - Q_ux'*inv_Q_uu*Q_ux;
         Vx(:,j)= Q_x - Q_ux'*inv_Q_uu*Q_u;
         V(:,j) = Q - 0.5*Q_u'*inv_Q_uu*Q_u;

    end 

    %----------------------------------------------> Find the controls
    % dx is initially zero because we start from the same point
    dx = zeros(2,1);    

    for i=1:(Horizon-1)    
         du = l_k(:,i) + L_k(:,:,i) * dx;   	% Feedback Controller 
         dx = A(:,:,i) * dx + B(:,:,i) * du;    % As we propagate forward, we use the linearized dynamics to approximate dx (this is the error from the nominal trajectory)
         u_new(:,i) = u_k(:,i) + gamma * du;    % Update controls with gamma to prevent controls from updating too fast
    end

    %Update nominal trajectory (u_k) for new updated controls
    u_k = u_new;    

    %---------------------------------------------> Simulation of the Nonlinear System
    %Create new nominal trajectory based on new control (u_new)
    [x_traj] = fnsimulate(xo,u_new,Horizon,dt,0);   
    [Cost(:,k)] =  fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);

    %fprintf('iLQG Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
 
end

u_traj = u_k;
