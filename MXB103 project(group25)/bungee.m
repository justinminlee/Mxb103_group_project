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

figure;

plot(0:dt:(T-dt), y);
xlabel('Time (s)');
ylabel('Distance from jumping point');
title('Position vs Time Using Modified Euler');

%% Plot: Velocity vs Time with Max Speed(Task 2)

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


%identifying index
D_val = H - D;
for i = 1:length(y) - 3
    if y(i) < D_val && y(i+1) < D_val && y(i+2) > D_val && y(i+3) > D_val
        break;
    end
end

% Time corresponding to each y value
t_vals = [(i-1)*dt, i*dt, (i+1)*dt, (i+2)*dt];
y_vals = [y(i), y(i+1), y(i+2), y(i+3)];

% Construct the Lagrange interpolating polynomial
p = @(t) sum(arrayfun(@(j) y_vals(j) * prod((t - t_vals([1:j-1, j+1:end])) ./ (t_vals(j) - t_vals([1:j-1, j+1:end]))), 1:4));


%Implement for bisection
a = t_vals(1);
b = t_vals(4);
tolerance = 1e-6;  %small value chosen for approximation.

while abs(b - a) > tolerance
    c = (a + b) / 2;
    if p(c) < D_val
        a = c;
    else
        b = c;
    end
end

t_camera_trigger = (a + b) / 2;

disp(['The camera should trigger at t = ', num2str(t_camera_trigger), ' seconds.']);


%% Task 6

% Initialization 
optimal_K = K;
optimal_L = L;
closest_distance = H - 1.75 - max(y);  % based on initial K, L

% Loop over a range of K values and L values to find the optimal parameters
for K_test = linspace(0.6*K, 1.4*K, 100) % Refine the search range
    for L_test = linspace(0.6*L, 1.4*L, 100) % Refine the search range
        [y_test, v_test, num_bounces_test, max_accel_test] = simulate_jump(K_test, L_test, dt, T, g, C, H);
        distance = H - 1.75 - max(y_test);
        
        % Prioritize water touch, then acceleration, and finally the number of bounces
        if distance >= 0 && abs(distance) < abs(closest_distance) && max_accel_test <= 2*g
            closest_distance = distance;
            optimal_K = K_test;
            optimal_L = L_test;
        end
    end
end

% Simulate jump with the optimal parameters
[y_optimal, v_optimal, num_bounces, max_accel_optimal] = simulate_jump(optimal_K, optimal_L, dt, T, g, C, H);

% Display the results
disp(['Closest distance to water: ', num2str(closest_distance), ' m']);
disp(['Optimal K(stiffness): ', num2str(optimal_K)]);
disp(['Optimal L(length): ', num2str(optimal_L)]);
disp(['Number of bounces: ', num2str(num_bounces)]);
disp(['Absolute maximum acceleration: ', num2str(max_accel_optimal), ' m/s^2']);
