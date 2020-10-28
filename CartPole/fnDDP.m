function  [u_traj, cost] = fnDDP(dynamicsf, x,u_new,num_iter,Horizon,gamma,p_target,dt,Q_f,R)

% Variables for the cart-pole problem
global g; 
global m_c;
global m_p; 
global l; 

g = 9.81;  		% gravity 
m_c = 1.0; 		% mass of the cart
m_p = 0.01;		% mass of the pole
l = 0.25;  		% length of the pole	

% Initialize the state space equations for the dynamics 
syms x1 x2 x3 x4 u

%Compute the Jacobians of the dynamics
dFx = matlabFunction(jacobian(dynamicsf, [x1; x2; x3; x4]));
dFu = matlabFunction(jacobian(dynamicsf, u));


% Initial State
xo = x;

% Initial Control:
u_k = u_new;
  
% Initial trajectory:
x_traj = zeros(4,Horizon);


    for k = 1:num_iter 	%Run for a certain number of iterations

        %------------------------------------------------> Linearization of the dynamics
        %------------------------------------------------> Quadratic Approximations of the cost function 
        for  j = 1:(Horizon-1) 	%Discretize trajectory for each timestep

            % linearization of dynamics (Jacobians dfx and dfu)
            [dfx,dfu] = Jacobians(x_traj(:,j),u_k(:,j),dFx, dFu);
            A(:,:,j) = eye(4,4) + dfx * dt; 	
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
        % initialize value function 
        Vxx(:,:,Horizon)= Q_f; 																		% Initialize Hessian of value function (Matrix)                                                                                         
        Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target);  										% Gradient of value function (Vector)
        V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target); 	% Value function (scalar)


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
        dx = zeros(4,1); 
        for i=1:(Horizon-1)    
           du = l_k(:,i) + L_k(:,:,i) * dx; 		%Feedback Controller 
           dx = A(:,:,i) * dx + B(:,:,i) * du; 		% As we propagate forward, we use the linearized dynamics to approximate dx (this is the error from the nominal trajectory)
           u_new(:,i) = u_k(:,i) + gamma * du; 		% Update controls with gamma to prevent controls from updating too fast
        end

        %Update nominal trajectory (u_k) for new updated controls
        u_k = u_new; 	


        %---------------------------------------------> Simulation of the Nonlinear System
        %Create new nominal trajectory based on new control (u_new)
        [x_traj] = fnsimulate(dynamicsf,xo,u_new,Horizon,dt,0);	 

    end

u_traj = u_k;
cost = fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);