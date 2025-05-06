clc; clear; close all;

load('F_external.mat');

% System parameters
m = 1.0;  % Mass (kg)
% Nonlinear spring parameters
k1 = 10.0;  % Linear coefficient (N/m)
k3 = 2.0;   % Cubic nonlinearity coefficient (N/m³)
% Nonlinear damper parameters
c1 = 0.5;   % Linear damping coefficient (N·s/m)
c3 = 0.05;  % Cubic nonlinearity coefficient (N·s³/m³)

% Nonlinear spring force function
F_spring = @(x) k1*x + k3*x.^3;

% Nonlinear damper force function
F_damper = @(v) c1*v + c3*v.^3;

function dxdt = system_dynamics(t, x, m, F_spring, F_damper, F_external)
    % Extract states
    x1 = x(1); x2 = x(2); x3 = x(3);
    v1 = x(4); v2 = x(5); v3 = x(6);
    
    % Calculate forces on each mass
    F1 = -F_spring(x1) - F_damper(v1) + F_spring(x2-x1) + F_damper(v2-v1);
    
    % Add external force to mass 2
    F2 = -F_spring(x2-x1) - F_damper(v2-v1) + F_spring(x3-x2) + F_damper(v3-v2) + F_external(t);
    
    F3 = -F_spring(x3-x2) - F_damper(v3-v2);
    
    % State derivatives
    dxdt = zeros(6,1);
    dxdt(1) = v1;
    dxdt(2) = v2;
    dxdt(3) = v3;
    dxdt(4) = F1/m;
    dxdt(5) = F2/m;
    dxdt(6) = F3/m;
end

% Load the external force functions
load('F_external.mat');

% Loop through each force type
for force_idx = 11:11 % length(F_external)
    % Get current force function
    current_force = F_external{force_idx};
    
    % Run simulation with current force
    disp(['Running simulation with force type ' num2str(force_idx)]);

    % Initial conditions
    x0 = [0, 0, 0, 0, 0, 0];  % Initial displacement of first mass only
    
    % Time span
    tspan = [0 300];  % Simulate for 60 seconds
    
    % ODE options
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    
    % Solve the system with current external force
    [t, x] = ode45(@(t,x) system_dynamics(t, x, m, F_spring, F_damper, current_force), tspan, x0, options);
end

% Plot positions
figure;
subplot(3,1,1);
plot(t, x(:,1), 'r-', t, x(:,2), 'g-', t, x(:,3), 'b-');
xlabel('Time (s)');
ylabel('Position (m)');
legend('Mass 1', 'Mass 2', 'Mass 3');
title('Mass Positions');
grid on;

% Plot velocities
subplot(3,1,2);
plot(t, x(:,4), 'r-', t, x(:,5), 'g-', t, x(:,6), 'b-');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('Mass 1', 'Mass 2', 'Mass 3');
title('Mass Velocities');
grid on;

% Plot external force
subplot(3,1,3);
% Evaluate the external force at each time point
F_ext_values = zeros(size(t));
for i = 1:length(t)
    F_ext_values(i) = F_external{force_idx}(t(i));
end
plot(t, F_ext_values, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Force (N)');
title('External Force on Mass 2');
grid on;