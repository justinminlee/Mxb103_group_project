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

% Integration using RK4 steps
for i = 1:N-1
   [y(i+1), v(i+1)] = RK4_steps(y(i), v(i), C, K, L, g, dt);
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

%% Plot: Position vs Time and velocity vs Time with Max Speed
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

%% Task 5 (look for 3.3 seconds for the answer(tutor mentioned)
% to be done

%% Task 6

% Step 1: Finding K and L for closest distance to water
optimal_K = K;
optimal_L = L;
closest_distance = H - 1.75 - min(y);  % based on initial K, L

for K_test = linspace(0.9*K, 1.1*K, 20)
    for L_test = linspace(0.9*L, 1.1*L, 20)
        [y_test, v_test] = simulate_jump(K_test, L_test, dt, T, g, C, H);
        distance = H - 1.75 - min(y_test);
        if distance < closest_distance
            closest_distance = distance;
            optimal_K = K_test;
            optimal_L = L_test;
        end
    end
end

% Step 2: Adjusting for number of bounces
desired_bounces = 10;
[y_optimal, v_optimal] = simulate_jump(optimal_K, optimal_L, dt, T, g, C, H);
sign_changes = find(diff(sign(v_optimal)));
num_bounces = length(sign_changes) / 2;

while abs(num_bounces - desired_bounces) > 1   % Allow 1 bounce difference for simplicity
    if num_bounces > desired_bounces
        optimal_L = 1.01 * optimal_L;  % Increase length slightly
    else
        optimal_L = 0.99 * optimal_L;  % Decrease length slightly
    end
    
    [y_optimal, v_optimal] = simulate_jump(optimal_K, optimal_L, dt, T, g, C, H);
    sign_changes = find(diff(sign(v_optimal)));
    num_bounces = length(sign_changes) / 2;
end

% Display the results
disp(['Optimal K: ', num2str(optimal_K)]);
disp(['Optimal L: ', num2str(optimal_L)]);
disp(['Number of bounces: ', num2str(num_bounces)]);


