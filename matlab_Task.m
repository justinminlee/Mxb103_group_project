clear; close all; clc;

% Parameters
H = 74; D = 31; c = 0.9; m = 80; L = 25; k = 90; g = 9.8;
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
    % K1
    K1_y = v(i);
    K1_v = g - C*abs(v(i))*v(i) - max(0, K*(y(i) - L));
    
    % K2
    K2_y = v(i) + 0.5*dt*K1_v;
    K2_v = g - C*abs(K2_y)*K2_y - max(0, K*(y(i) + 0.5*dt*K1_y - L));
    
    % K3
    K3_y = v(i) + 0.5*dt*K2_v;
    K3_v = g - C*abs(K3_y)*K3_y - max(0, K*(y(i) + 0.5*dt*K2_y - L));
    
    % K4
    K4_y = v(i) + dt*K3_v;
    K4_v = g - C*abs(K4_y)*K4_y - max(0, K*(y(i) + dt*K3_y - L));
    
    % Update
    y(i+1) = y(i) + (dt/6)*(K1_y + 2*K2_y + 2*K3_y + K4_y);
    v(i+1) = v(i) + (dt/6)*(K1_v + 2*K2_v + 2*K3_v + K4_v);
end




% Analysis: Number of Bounces and Timing
[~,locs] = findpeaks(-y); 
bounces = length(locs);
time_for_bounces = locs(end) * dt;

fprintf('Number of bounces: %d\n', bounces);
fprintf('Time for bounces: %.2f s\n', time_for_bounces);

% Analysis: Maximum Speed and Timing
[max_speed, max_speed_idx] = min(v);
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

plot(0:dt:(T-dt), v);
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
