% After running your simulation and obtaining t and x
figure;

% Mass 1 Phase Portrait with Vector Field
subplot(1,3,1);
[X1, V1] = meshgrid(linspace(min(x(:,1))*1.2, max(x(:,1))*1.2, 20), ...
                    linspace(min(x(:,4))*1.2, max(x(:,4))*1.2, 20));
DX1 = V1;  % dx/dt = v
DV1 = zeros(size(X1));

% Calculate acceleration field for a mesh of points for Mass 1
for i = 1:size(X1, 1)
    for j = 1:size(X1, 2)
        % Simplified calculation for visualization only
        spring_force = F_spring(X1(i,j));
        damper_force = F_damper(V1(i,j));
        DV1(i,j) = (-spring_force - damper_force)/m;
    end
end

% Plot vector field and trajectory for Mass 1
quiver(X1, V1, DX1, DV1, 0.5, 'b');
hold on;
plot(x(:,1), x(:,4), 'r-', 'LineWidth', 1);
% Mark the first point with a bold black dot
plot(x(1,1), x(1,4), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% Mark the last point with a thin black dot
plot(x(end,1), x(end,4), 'k*', 'MarkerSize', 10, 'MarkerFaceColor', 'none');
xlabel('Position (m)');
ylabel('Velocity (m/s)');
title('Mass 1 Phase Portrait');
grid on;
axis square;

% Mass 2 Phase Portrait with Vector Field
subplot(1,3,2);
[X2, V2] = meshgrid(linspace(min(x(:,2))*1.2, max(x(:,2))*1.2, 20), ...
                    linspace(min(x(:,5))*1.2, max(x(:,5))*1.2, 20));
DX2 = V2;  % dx/dt = v
DV2 = zeros(size(X2));

% Calculate acceleration field for a mesh of points for Mass 2
for i = 1:size(X2, 1)
    for j = 1:size(X2, 2)
        % Simplified calculation for visualization only
        spring_force = F_spring(X2(i,j));
        damper_force = F_damper(V2(i,j));
        DV2(i,j) = (-spring_force - damper_force)/m;
    end
end

% Plot vector field and trajectory for Mass 2
quiver(X2, V2, DX2, DV2, 0.5, 'b');
hold on;
plot(x(:,2), x(:,5), 'r-', 'LineWidth', 1);
% Mark the first point with a bold black dot
plot(x(1,2), x(1,5), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% Mark the last point with a thin black dot
plot(x(end,2), x(end,5), 'k*', 'MarkerSize', 10);
xlabel('Position (m)');
ylabel('Velocity (m/s)');
title('Mass 2 Phase Portrait');
grid on;
axis square;

% Mass 3 Phase Portrait with Vector Field
subplot(1,3,3);
[X3, V3] = meshgrid(linspace(min(x(:,3))*1.2, max(x(:,3))*1.2, 20), ...
                    linspace(min(x(:,6))*1.2, max(x(:,6))*1.2, 20));
DX3 = V3;  % dx/dt = v
DV3 = zeros(size(X3));

% Calculate acceleration field for a mesh of points for Mass 3
for i = 1:size(X3, 1)
    for j = 1:size(X3, 2)
        % Simplified calculation for visualization only
        spring_force = F_spring(X3(i,j));
        damper_force = F_damper(V3(i,j));
        DV3(i,j) = (-spring_force - damper_force)/m;
    end
end

% Plot vector field and trajectory for Mass 3
quiver(X3, V3, DX3, DV3, 0.5, 'b');
hold on;
plot(x(:,3), x(:,6), 'r-', 'LineWidth', 1);
% Mark the first point with a bold black dot
plot(x(1,3), x(1,6), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% Mark the last point with a thin black dot
plot(x(end,3), x(end,6), 'k*', 'MarkerSize', 10);
xlabel('Position (m)');
ylabel('Velocity (m/s)');
title('Mass 3 Phase Portrait');
grid on;
axis square;

% Adjust the figure to look better
set(gcf, 'Position', [100, 100, 1200, 400]);  % Width is 3x height for better horizontal display