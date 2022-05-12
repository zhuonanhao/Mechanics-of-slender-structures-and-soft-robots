%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Filename: Elastic_Beam_Bending.m
%   Description: Simulation of the deformation of elastic beams and 
%                comparison with Euler-Bernoulli beam theory
%   Author: Zhuonan Hao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
close all;

%% Parameters 

%-----------------------------------------------------------
%   Simulation parameters
%-----------------------------------------------------------

N = 50; % Number of vertices
ne = N-1; % Number of edges
dt = 0.01; % Timestep size, second
totalTime = 1; % Total time, second
Nsteps = round(totalTime/dt); % Number of steps

%-----------------------------------------------------------
%   Geometric parameters
%-----------------------------------------------------------

R = 0.013; % Outer radius, m
r = 0.011; % Inner radius, m
A = pi*(R^2-r^2); % Cross section area, m^2

l = 1; % Rod length, m
deltaL = l/(N-1); % Discrete length, m

% Structure definition
nodes = zeros(N,2);
for kk = 1:N
    nodes(kk,1) = (kk-1)*deltaL;
    nodes(kk,2) = 0;
end

free_index = 3:2*N-1; % Free index
fixed_index = [1,2,2*N]; % Fixed index

%-----------------------------------------------------------
%   Material parameters
%-----------------------------------------------------------
Y = 70e9; % Young's modulus, N/m^2(Pa)
I = pi/4*(R^4-r^4); % Moment of inertia of the cross section, m^4

rho = 2700; % Density, kg/m^3

mi = pi*(R^2-r^2)*l*rho/(N-1); % Node Mass, kg

% Mass matrix
M = zeros(2*N,2*N);
for kk = 1:N
    M(2*kk-1,2*kk-1) = mi;
    M(2*kk,2*kk) = mi;
end

EI = Y*I; % Bending stiffness, N-m^2
EA = Y*A; % Stretching stiffness, N

%-----------------------------------------------------------
%   Physical parameters
%-----------------------------------------------------------
g = 9.8; % Gravitational acceleration, N/m^2

P = 2000; % Applied force, N
d = 0.75; % Distance of the applied force to the left end, m 

% Find the index of the node with applied force
[~,nodes_force_index] = min(abs(nodes(:,1)-d)); 

% External force, N
W = zeros(2*N,1);
for kk = 1:N
    W(2*kk) = -mi*g;
    if kk == nodes_force_index
        W(2*kk) = W(2*kk) - P;
    end
end

%% Initialization

% DOF Vecotr
q0 = zeros(2*N,1);
for kk = 1:N
    q0(2*kk-1) = nodes(kk,1); % x coordinate 
    q0(2*kk) = nodes(kk,2); % y coordinate
end

% Initial condition
q = q0; % Position, m
u = (q-q0)/dt; % Velocity, m/s

% Solutions
q_solution = zeros(Nsteps,2*N); % Position, m
u_solution = zeros(Nsteps,2*N); % Velocity, m/s
q_solution(1,:) = q0';
u_solution(1,:) = (q_solution(1,:)-q0')/dt; 

%% Simulation

for kk = 2:Nsteps
    
    % Time marching
    time = (kk-1)*dt;
    fprintf('Time = %f\n', time);
    
    
    %                   Newton-Raphson method
    %-----------------------------------------------------------
    %   Define tolerence and error
    %-----------------------------------------------------------
    tol = EI/l^2*1e-3; % Tolerence
    err = 10*tol; % Error
    
    %-----------------------------------------------------------
    %   Guess solution
    %-----------------------------------------------------------
    q = q0;
    
    %-----------------------------------------------------------
    %   Update the fixed index
    %-----------------------------------------------------------
    q_fixed = [0;0;0];
    q(fixed_index) = q_fixed;
    
    %-----------------------------------------------------------
    %   Update the free index
    %-----------------------------------------------------------    
    q_free = q(free_index);
    
    while err > tol

        %-----------------------------------------------------------
        %   Compute F and J
        %-----------------------------------------------------------

        % Inertia forces
        f = M/dt*((q-q0)/dt-u);
        J = M/dt^2;

        % Elastic forces - streching force between nodes k and k+1
        for k = 1:N-1

            % Index
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);

            % Compute 
            dF = gradEs(xk,yk,xkp1,ykp1,deltaL,EA);
            dJ = hessEs(xk,yk,xkp1,ykp1,deltaL,EA);

            % Add
            f(2*k-1:2*k+2) = f(2*k-1:2*k+2)+dF;
            J(2*k-1:2*k+2,2*k-1:2*k+2) = J(2*k-1:2*k+2,2*k-1:2*k+2)+dJ;
        end

        % Elastic forces - bending force between nodes k-1, k, and k+1
        for k = 2:N-1

            % Index
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);

            curvature0 = 0;

            % Compute 
            dF = gradEb(xkm1,ykm1,xk,yk,xkp1,ykp1,curvature0,deltaL,EI);
            dJ = hessEb(xkm1,ykm1,xk,yk,xkp1,ykp1,curvature0,deltaL,EI);

            % Add
            f(2*k-3:2*k+2) = f(2*k-3:2*k+2)+dF;
            J(2*k-3:2*k+2,2*k-3:2*k+2) = J(2*k-3:2*k+2,2*k-3:2*k+2)+dJ;
        end

        % External forces 
        f = f-W;

        %-----------------------------------------------------------
        %   Update q
        %-----------------------------------------------------------
        f_free = f(free_index);
        J_free = J(free_index,free_index);
        q_free = q_free-J_free\f_free;
        
        %-----------------------------------------------------------
        %   Compute error
        %-----------------------------------------------------------
        err = sum(abs(f_free));
        
        %-----------------------------------------------------------
        % Plug free DOFs back into the full DOF vector
        %-----------------------------------------------------------
        q(free_index) = q_free;
    end
    
    % Update position and volocity
    u = (q-q0)/dt;
    q0 = q;
    
    % Store solution
    q_solution(kk,:) = q';
    u_solution(kk,:) = u';
    
    % Animation
    if mod(time,dt) == 0
        figure(1)
        plot(q(1:2:end),q(2:2:end),'ro-')
        axis equal 
        drawnow
    end
end
