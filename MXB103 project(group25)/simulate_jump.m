function [y, v, num_bounces, max_acceleration] = simulate_jump(K, L, dt, T, g, C, H)
    % Initialization
    N = T/dt;
    y = zeros(1, N);
    v = zeros(1, N);
    y(1) = 0; % Initial position
    v(1) = 0; % Initial velocity

    max_acceleration = 0;  % to store the maximum acceleration value during the jump

    % RK4 calculations using your function
    for i = 1:N-1
       [y(i+1), v(i+1)] = RK4_steps(y(i), v(i), C, K, L, g, dt);
       acceleration = g - C*abs(v(i+1))*v(i+1) - max(0,K*(y(i+1)-L));  % compute acceleration at each step
       max_acceleration = max(max_acceleration, abs(acceleration));  % update if it's a new max
    end

    % Calculate bounces
    sign_changes = find(diff(sign(v)));
    num_bounces = length(sign_changes) / 2;
end

