%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Filename: Three_Node_Implicit.m
%   Description: Simulation of the motion of a sphere falling inside 
%                viscous fluid
%   Author: Zhuonan Hao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
close all;

%% Parameters 

%-----------------------------------------------------------
%   Simulation parameters
%-----------------------------------------------------------

N = 3; % Number of vertices
ne = N-1; % Number of edges
dt = 0.01; % Timestep size, second
totalTime = 10; % Total time, second
Nsteps = round(totalTime/dt); % Number of steps

%-----------------------------------------------------------
%   Geometric parameters
%-----------------------------------------------------------
% Radius of spheres, m
R1 = 0.005;
R2 = 0.025;
R3 = 0.005;

r0 = 0.001; % Rod radius, m
RodLength = 0.1; % Rod length, m
deltaL = RodLength/(N-1); % Discrete length, m

% Structure definition
nodes = zeros(N,2);
for kk = 1:N
    nodes(kk,1) = (kk-1)*deltaL;
    nodes(kk,2) = 0;
end

%-----------------------------------------------------------
%   Material parameters
%-----------------------------------------------------------

visc = 1000; % Viscosity, Pa-s
Y = 1e9; % Young's modulus, N/m^2(Pa)

% Density, kg/m^3
rho_metal = 7000;
rho_fluid = 1000;

% Viscous damping, N/m
C1 = 6*pi*visc*R1;
C2 = 6*pi*visc*R2;
C3 = 6*pi*visc*R3;

% Viscous damping matrix
C = diag([C1,C1,C2,C2,C3,C3]);

% Mass, kg
m1 = 4/3*pi*R1^3*rho_metal;
m2 = 4/3*pi*R2^3*rho_metal;
m3 = 4/3*pi*R3^3*rho_metal;

% Mass matrix
M = diag([m1,m1,m2,m2,m3,m3]);

EI = Y*pi*r0^4/4; % Bending stiffness, N-m^2
EA = Y*pi*r0^2; % Stretching stiffness, N

%-----------------------------------------------------------
%   Physical parameters
%-----------------------------------------------------------
g = 9.8; % Gravitational acceleration, N/m^2

% External force, N
W = zeros(2*N,1);
W(2) = -4/3*pi*R1^3*(rho_metal-rho_fluid)*g;
W(4) = -4/3*pi*R2^3*(rho_metal-rho_fluid)*g;
W(6) = -4/3*pi*R3^3*(rho_metal-rho_fluid)*g;

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
    tol = EI/RodLength^2*1e-3; % Tolerence
    err = 10*tol; % Error
    
    %-----------------------------------------------------------
    %   Guess solution
    %-----------------------------------------------------------
    q = q0;
    
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

        % Viscous forces
        f = f+C*(q-q0)/dt;
        J = J+C/dt;

        % External forces 
        f = f-W;

        %-----------------------------------------------------------
        %   Update q
        %-----------------------------------------------------------
        q = q-J\f;
        
        %-----------------------------------------------------------
        %   Compute error
        %-----------------------------------------------------------
        err = sum(abs(f));
        
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



 














