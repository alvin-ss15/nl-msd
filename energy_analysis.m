% Calculate energies at each time step
kinetic_energy = zeros(length(t), 3);  % KE for each mass
potential_energy = zeros(length(t), 3);  % PE for each spring
damper_energy = zeros(length(t), 3);  % Energy dissipated by each damper
total_energy = zeros(length(t), 1);    % Total system energy

% Calculate spring potential energy integral function
% For nonlinear spring F(x) = k1*x + k3*x^3
% PE = ∫F(x)dx = (1/2)*k1*x^2 + (1/4)*k3*x^4
PE_spring = @(x) (1/2)*k1*x.^2 + (1/4)*k3*x.^4;

for i = 1:length(t)
    % Positions and velocities at time step i
    x1 = x(i,1); x2 = x(i,2); x3 = x(i,3);
    v1 = x(i,4); v2 = x(i,5); v3 = x(i,6);
    
    % Kinetic energy: KE = (1/2)*m*v^2
    kinetic_energy(i,1) = 0.5*m*v1^2;
    kinetic_energy(i,2) = 0.5*m*v2^2;
    kinetic_energy(i,3) = 0.5*m*v3^2;
    
    % Potential energy in springs
    potential_energy(i,1) = PE_spring(x1);  % Wall to mass 1
    potential_energy(i,2) = PE_spring(x2-x1);  % Mass 1 to mass 2
    potential_energy(i,3) = PE_spring(x3-x2);  % Mass 2 to mass 3
    
    % Calculate energy dissipated by dampers
    % This requires integration over time, here's a simplified approach
    % using trapezoidal rule if i > 1
    if i > 1
        dt = t(i) - t(i-1);
        
        % Damper velocities
        v_d1 = v1;  % Wall to mass 1
        v_d2 = v2 - v1;  % Mass 1 to mass 2
        v_d3 = v3 - v2;  % Mass 2 to mass 3
        
        % Damper forces
        f_d1 = F_damper(v_d1);
        f_d2 = F_damper(v_d2);
        f_d3 = F_damper(v_d3);
        
        % Power dissipated = Force * Velocity
        damper_energy(i,1) = damper_energy(i-1,1) + f_d1 * v_d1 * dt;
        damper_energy(i,2) = damper_energy(i-1,2) + f_d2 * v_d2 * dt;
        damper_energy(i,3) = damper_energy(i-1,3) + f_d3 * v_d3 * dt;
    end
    
    % Total energy = KE + PE (excluding dissipated energy)
    total_energy(i) = sum(kinetic_energy(i,:)) + sum(potential_energy(i,:));
end

% Plot energy components
figure;

% Kinetic energy
subplot(4,1,1);
plot(t, kinetic_energy);
xlabel('Time (s)');
ylabel('Kinetic Energy (J)');
legend('Mass 1', 'Mass 2', 'Mass 3');
title('Kinetic Energy');
grid on;

% Potential energy
subplot(4,1,2);
plot(t, potential_energy);
xlabel('Time (s)');
ylabel('Potential Energy (J)');
legend('Spring 1', 'Spring 2', 'Spring 3');
title('Potential Energy');
grid on;

% Damper dissipated energy
subplot(4,1,3);
plot(t, damper_energy);
xlabel('Time (s)');
ylabel('Dissipated Energy (J)');
legend('Damper 1', 'Damper 2', 'Damper 3');
title('Energy Dissipated by Dampers');
grid on;

% Total mechanical energy
subplot(4,1,4);
plot(t, total_energy, 'k-', ...
     t, total_energy + sum(damper_energy, 2), 'r--');
xlabel('Time (s)');
ylabel('Energy (J)');
legend('Mechanical Energy', 'Total Energy (incl. dissipated)');
title('System Energy');
grid on;