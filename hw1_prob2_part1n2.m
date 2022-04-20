%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Filename: N_Node_Implicit.m
%   Description: Simulation of the motion of N-connected spheres falling 
%                inside viscous fluid
%   Author: Zhuonan Hao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
close all;

%% Parameters 

%-----------------------------------------------------------
%   Simulation parameters
%-----------------------------------------------------------

N = 21; % Number of vertices
ne = N-1; % Number of edges
dt = 0.01; % Timestep size, second
totalTime = 50; % Total time, second
Nsteps = round(totalTime/dt); % Number of steps

%-----------------------------------------------------------
%   Geometric parameters
%-----------------------------------------------------------

r0 = 0.001; % Rod radius, m
RodLength = 0.1; % Rod length, m
deltaL = RodLength/(N-1); % Discrete length, m

% Radius of spheres, m
R = ones(N,1)*deltaL/10;
middle_index = ceil(N/2);
R(middle_index) = 0.025;

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
Ci = 6*pi*visc*R;

% Viscous damping matrix
C = zeros(2*N,2*N);
for kk = 1:N
    C(2*kk-1,2*kk-1) = Ci(kk);
    C(2*kk,2*kk) = Ci(kk);
end

% Mass, kg
mi = 4/3*pi*R.^3*rho_metal;

% Mass matrix
M = zeros(2*N,2*N);
for kk = 1:N
    M(2*kk-1,2*kk-1) = mi(kk);
    M(2*kk,2*kk) = mi(kk);
end

EI = Y*pi*r0^4/4; % Bending stiffness, N-m^2
EA = Y*pi*r0^2; % Stretching stiffness, N

%-----------------------------------------------------------
%   Physical parameters
%-----------------------------------------------------------
g = 9.8; % Gravitational acceleration, N/m^2

% External force, N
W = zeros(2*N,1);
for kk = 1:N
    W(2*kk) = -4/3*pi*R(kk)^3*(rho_metal-rho_fluid)*g;
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
    
end

%% Results analysis

%-----------------------------------------------------------
%   Position and velocity of the middle node
%-----------------------------------------------------------

middle_position_y = q_solution(:,middle_index*2);
middle_velocity_y = u_solution(:,middle_index*2);

time = 0:dt:totalTime-dt;

figure(1)
plot(time,middle_position_y,'linewidth',3)
title('Middle node')
xlabel('Time, t [s]')
ylabel('Position, y [m]')
set(gca,'fontsize',16)

figure(2)
plot(time,middle_velocity_y,'linewidth',3)
title('Middle node')
xlabel('Time, t [s]')
ylabel('Velocity, dy/dt [m/s]')
set(gca,'fontsize',16)
 
%-----------------------------------------------------------
%   Final deformed shape
%-----------------------------------------------------------

x_final = q_solution(end,1:2:end);
y_final = q_solution(end,2:2:end);

figure(3)
plot(x_final,y_final,'ro-','linewidth',3)
title('Final deformed shape')
xlabel('Position, x [m]')
ylabel('Position, y [m]')
set(gca,'fontsize',16)













