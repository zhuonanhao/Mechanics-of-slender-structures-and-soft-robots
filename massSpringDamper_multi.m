% Filename: massSpringDamper_multi.m
% Simulation of multi-mass, spring, damper system

clear all;
close all;

N = 3; % number of DOF

m = [10; 0.01; 0.001]; %#ok<*NASGU> % mass
K = [1; 5; 0.1]; % spring constant
b = [10; 50; 100]; % viscous damping coefficient
F0 = [0.1; 10; 1]; % driving amplitude
omega = [0.1; 10; 2]; % driving frequency
x0 = [1; 2; 3];

maxTime = 10; % total time of simulation
dt = 1e-2; % discrete time step
eps = 1e-6 * sum(F0(:)); % error tolerance

t = 0:dt:maxTime; % time
x = zeros(N, numel(t)); % position
u = zeros(N, numel(t)); % velocity

% Initial condition: position and velocity are zero
x(:, 1) = 0;
u(:, 1) = 0;

f = zeros(N,1);
J = zeros(N,N);

for k=1:length(t)-1 % march over time steps
    x_old = x(:, k);
    u_old = u(:, k);    
    t_new = t(k+1);
    
    x_new = x_old; % guess solution
    err = eps * 100; % initialize to a large value
    while err > eps
        % Compute f
        computeF
        % Compute J
        computeJ
        % Compute delta x
        deltaX = J \ f;
        x_new = x_new - deltaX;
        err = abs( f );
    end
    
    u_new = (x_new - x_old) / dt;
    
    x(:, k+1) = x_new; % store solution
    u(:, k+1) = u_new;
    
    % Plot
    figure(1);
    plot([0; x_new + x0], zeros(N+1,1), 'ro-');
    xlim([0 max(x0)*1.5]);
    ylim([-1 1]);
    title(num2str(t_new, 't=%5.2f sec'));
    axis off
    drawnow
end

% Plot it
FONT = 'Arial';
FONTSIZE = 10;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81;237 28 36;0 174 239; 0 0 0]/255; % colors

h1 = figure(1);
plot(t, x(1,:), 'Color', colpos(1,:), 'LineWidth', 2);
hold on
plot(t, x(2,:), 'Color', colpos(2,:), 'LineWidth', 2);
plot(t, x(3,:), 'Color', colpos(3,:), 'LineWidth', 2);
hold off
box on
xlabel('Time, t [sec]','Fontname', FONT,'FontSize',FONTSIZE);
ylabel('Location, x_i [m]','Fontname', FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
l = legend('i=1','i=2','i=3');
set(l,'Fontname', FONT,'FontSize',FONTSIZE);
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
saveas(h1, 'Fig_multiMassSpringDamper.pdf');
