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
y(1) = 0; % Initial position

%% Classical Fourth Order Runge-Kutta (RK4) Method
%Task 1
% definitely has to be a function
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





  %  Adjust the initial conditions to start from an equilibrium position.
  %  Use an initial velocity to simulate the upward motion.
  %  Adjust damping to ensure the amplitude decreases over time.



% Analysis: Number of Bounces and Timing
[~,locs] = findpeaks(-y); 
bounces = length(locs);
time_for_bounces = locs(end) * dt;

fprintf('Number of bounces: %d\n', bounces);
fprintf('Time for bounces: %.2f s\n', time_for_bounces);

% Analysis: Maximum Speed and Timing
[max_speed, max_speed_idx] = max(v); % we are supposed to find the max, so max(v) not min(v)
time_of_max_speed = max_speed_idx * dt;

fprintf('Maximum speed: %.2f m/s\n', abs(max_speed));
fprintf('Time of maximum speed: %.2f s\n', time_of_max_speed);

%% Plot: Position vs Time(Task 1)
%this is correct(looks correct)
figure;

plot(0:dt:(T-dt), y);
xlabel('Time (s)');
ylabel('Distance from jumping point');
title('Position vs Time Using Modified Euler');

%% Plot: Velocity vs Time with Max Speed(Task 2)
%(fixed this as well)
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

%% Task 3
% Initialize acceleration array
a = zeros(1, N);

% Calculate acceleration using central difference for all points except the first and last
for i = 2:N-1
    a(i) = (v(i+1) - v(i-1)) / (2 * dt);
end

% For the first and last points, use forward and backward difference
a(1) = (v(2) - v(1)) / dt;
a(N) = (v(N) - v(N-1)) / dt;


% Find the max positive and max negative acceleration and their timings
[max_pos_accel, max_pos_accel_idx] = max(a);
[max_neg_accel, max_neg_accel_idx] = min(a);

% Determine which acceleration is of greater magnitude
if abs(max_pos_accel) > abs(max_neg_accel)
    max_accel = max_pos_accel;
    max_accel_idx = max_pos_accel_idx;
else
    max_accel = max_neg_accel;
    max_accel_idx = max_neg_accel_idx;
end
time_of_max_accel = max_accel_idx * dt;

% Plot acceleration vs time
figure;
plot(0:dt:(T-dt), a);
hold on;
plot(time_of_max_accel, max_accel, 'ro');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
title('Acceleration vs Time');
legend('Acceleration', 'Max Acceleration');
hold off;

% Print maximum acceleration and its timing
fprintf('Maximum acceleration: %.2f m/s^2\n', abs(max_accel));
fprintf('Time of maximum acceleration: %.2f s\n', time_of_max_accel);

%% Task 4

% Calculate the integral using the trapezoidal rule
total_distance = sum(abs(v) + [0, abs(v(1:end-1))]) * 0.5 * dt;

fprintf('Total distance traveled by the jumper in 60 seconds: %.2f m\n', total_distance);

%% Task 5 (look for 3.3 seconds for the answers)

% 1. Find the Nearest Four Values
index = find(y < (H-D), 2, 'last');  % Last two values below H-D
index = [index; find(y > (H-D), 2, 'first')]; % First two values above H-D

% 2. Interpolation
t_values = index*dt;  % Convert indices to time values
y_values = y(index);

% Constructing the polynomial using MATLAB's polyfit
p_coefficients = polyfit(t_values, y_values, 3); % 3rd degree polynomial

% 3. Root Finding
% Create a function handle for the polynomial minus (H-D)
p = @(t) polyval(p_coefficients, t) - (H-D);

% Using MATLAB's fzero for root finding
time_for_camera = fzero(p, [t_values(2), t_values(3)]); 

fprintf('The camera should trigger at: %.2f s\n', time_for_camera);

