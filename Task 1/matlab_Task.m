
clear; close all; clc;

% Parameters
H = 74; % Height of jump point
D = 31; % Deck height
c = 0.9; % Drag coefficient 
m = 80; % Mass of jumper
L = 25; % Length of bungee rope
k = 90; % Spring constant of bungee rope
g = 9.8; % Gravitational acceleration 

%scaled parameters

C = c/m;
K = k/m;
dt = 0.01; % time step
T = 60; % total time
N = T/dt; % number of steps
y = zeros(1, N); v = zeros(1, N); % Initialize position and velocity
y(1) = H; % Initial position

%% Classical Fourth Order Runge-Kutta (RK4) Method
for i = 1:N-1
    [y(i+1), v(i+1)] = RK4_step(y(i), v(i), C, K, L, g, dt);
end




% Analysis: Number of Bounces and Timing
[~,locs] = findpeaks(-y); 
bounces = length(locs);
time_for_bounces = locs(end) * dt;

fprintf('Number of bounces: %d\n', bounces);
fprintf('Time for bounces: %.2f s\n', time_for_bounces);

% Analysis: Maximum Speed and Timing
[max_speed, max_speed_idx] = max(v);
time_of_max_speed = max_speed_idx * dt;

fprintf('Maximum speed: %.2f m/s\n', abs(max_speed));
fprintf('Time of maximum speed: %.2f s\n', time_of_max_speed);

%% Plot: Position vs Time
figure;

plot(0:dt:(T-dt), y);
xlabel('Time (s)');
ylabel('Position (m)');
title('Position vs Time Using Modified Euler');

%% Plot: Velocity vs Time with Max Speed
figure;

plot(0:dt:(T-dt), -v);
hold on;
plot(time_of_max_speed, max_speed, 'ro');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs Time');
legend('Velocity', 'Max Speed');
hold off;

%% Plot: Position vs Time and elocity vs Time with Max Speed
figure;
plot(0:dt:(T-dt), y);
xlabel('Time (s)');
ylabel('Position (m)-Velocity (m/s)');
title('Position/Velocity vs Time Using Modified Euler');
hold on;
plot(0:dt:(T-dt), v);

plot(time_of_max_speed, max_speed, 'ro');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs Time');
legend('Pos','Velocity', 'Max Speed');
hold off;
