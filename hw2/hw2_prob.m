%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Filename: DiscreteElasticRods.m
%   Description: Simulation of the elastic rod
%   Author: Zhuonan Hao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
close all;

%% Parameters 

%-----------------------------------------------------------
%   Simulation parameters
%-----------------------------------------------------------

nv = 50; % Number of vertices
ne = nv-1; % Number of edges
dt = 0.01; % Timestep size, second
totalTime = 5; % Total time, second
maximum_iter = 100; % Maximum iterations in Newton's solver
tol = 1e-3; % Tolerance on force function (to be multipled by characteristic bending force)
ndof = 4*nv-1; % Number of degrees of freedom
Nsteps = round(totalTime/dt); % Number of steps

%-----------------------------------------------------------
%   Geometric parameters
%-----------------------------------------------------------

r0 = 0.001; % Rod cross-section radius, m
RodLength = 0.2; % Rod length, m
natR = 0.02; % Natural radius of curvature

% Structure definition
nodes = zeros(nv,3);
dTheta = (RodLength/ne)/natR;
for kk = 1:nv
    nodes(kk,1) = natR*cos((kk-1)*dTheta);
    nodes(kk,2) = natR*sin((kk-1)*dTheta);
end

% Reference length
refLen = zeros(ne,1);
for kk = 1:ne
    dx = nodes(kk+1,:)-nodes(kk,:);
    refLen(kk) = norm(dx);
end

% Voronoi length
voronoiRefLen = zeros(nv,1);
for kk = 1:nv
    if kk == 1
        voronoiRefLen(kk) = 0.5*refLen(kk);
    elseif kk == nv
        voronoiRefLen(kk) = 0.5*refLen(kk-1);
    else
        voronoiRefLen(kk) = 0.5*(refLen(kk-1)+refLen(kk));
    end
end

% Reference directors
d1 = zeros(ne,3);
d2 = zeros(ne,3);

tangent = zeros(ne,3);
for kk = 1:ne
    dx = nodes(kk+1,:)-nodes(kk,:);
    tangent(kk,:) = dx/norm(dx);
end

% Figure out a good choice for d1(1,:)
t0 = tangent(1,:);
t1 = [0 0 -1];
d1Tmp = cross(t0,t1);
if abs(d1Tmp) < 1e-6
    t1 = [0 1 0];
    d1Tmp = cross(t0,t1);
end
d1(1,:) = d1Tmp/norm(d1Tmp);
d2Tmp = cross(t0,d1(1,:));
d2(1,:) = d2Tmp/norm(d2Tmp);

% Parallel transport for edge 2 to edge ne
for kk = 2:ne
    t0 = tangent(kk-1,:);
    t1 = tangent(kk,:);
    d1_old = d1(kk-1,:);
    d1_new = parallel_transport(d1_old,t0,t1);
    d1(kk,:) = d1_new/norm(d1_new);
    d2_new = cross(t1,d1_new);
    d2(kk,:) = d2_new/norm(d2_new);
end

%-----------------------------------------------------------
%   Material parameters
%-----------------------------------------------------------

rho = 1000; % Density, kg/m^3
nu = 0.5; % Poisson ratio
Y = 1e9; % Young's modulus, N/m^2(Pa)
G = Y/(2*(1+nu)); % Shear modulus, N/m^2(Pa)
EI = Y*pi*r0^4/4; % Bending stiffness, N-m^2
EA = Y*pi*r0^2; % Stretching stiffness, N
GJ = G*pi*r0^4/2; % Twisting stiffness, N-m^2
dm = (pi*r0^2*RodLength)*rho/ne; % Unit mass, kg

% Mass matrix
m = zeros(ndof,1);
for kk = 1:nv
    if kk == 1 || kk == nv
        m(4*(kk-1)+1:4*(kk-1)+3) = dm/2;
    else
        m(4*(kk-1)+1:4*(kk-1)+3) = dm;
    end
end
for kk = 1:ne
    m(4*kk) = dm/2*r0^2;
end

%-----------------------------------------------------------
%   Physical parameters
%-----------------------------------------------------------

g = [0;0;-9.81]; % Gravitational acceleration, N/m^2
garr = zeros(ndof,1);
for kk = 1:nv
    garr(4*(kk-1)+1:4*(kk-1)+3) = g;
end
Fg = m.*garr; % Gravity vector

%% Initialization

% DOF vector
q0 = zeros(ndof,1); 
for kk = 1:nv
    q0(4*(kk-1)+1:4*(kk-1)+3) = nodes(kk,:); % x,y,z
end
q0(4:4:end) = 0; % theta

% Initial condition
q = q0;
u = (q-q0)/dt;

% Constrained and unconstrained indices
consInd = 1:7; % Fixed 
UnconsInd = 8:ndof; % Free

% Material director
theta = q0(4:4:end);
[m1,m2] = computeMaterialDirectors(d1,d2,theta);

% Reference twist
refTwist = zeros(nv,1);
refTwist = getRefTwist(d1,tangent,refTwist);

% Natural curvature
kappaBar = getkappa(q,m1,m2);

%% Simulation
endZ = zeros(Nsteps, 1); % z-coordinate of the last node
currentTime = 0; % Current time

d1_old = d1;
d2_old = d2;
refTwist_old = refTwist;

figure(1)
hold on
for timeStep = 1:Nsteps 
    
    [q, d1, d2, refTwist] = DER(q0, u, d1_old, d2_old, refTwist_old); % Main function that steps forward one time step

    % Plot
    x_coord = q(1:4:end);
    y_coord = q(2:4:end);
    z_coord = q(3:4:end);
    figure(1);
    clf();
    plot3(x_coord, y_coord, z_coord, 'ro-');
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    xlim([-0.05,0.05]);
    ylim([-0.05,0.05]);
    zlim([-0.1,0.05]);
    box on
    endZ(timeStep) = q(end);

    u = (q - q0) / dt; % Update velocity
    q0 = q; % New position becomes old position
    d1_old = d1; % New reference director becomes old reference director
    d2_old = d2;
    refTwist_old = refTwist; % New reference twist becomes old reference twist
    currentTime = currentTime + dt; % Current time
end

%% Result Analysis

%-----------------------------------------------------------
%   z-coordinate of the last node
%-----------------------------------------------------------
figure(2)
tarray = (1:Nsteps) * dt;
plot(tarray, endZ, 'r-', 'linewidth',3);
xlabel('Time, t [s]');
ylabel('Tip displacement, \delta_z [m]');
set(gca,'fontsize',16)



