function [y, v] = simulate_jump(K, L, dt, T, g, C, H)
    % Initialization
    N = T/dt;
    y = zeros(1, N);
    v = zeros(1, N);
    y(1) = 0; % Initial position
    v(1) = 0; % Initial velocity

    % RK4 calculations using your function
    for i = 1:N-1
       [y(i+1), v(i+1)] = RK4_steps(y(i), v(i), C, K, L, g, dt);
    end
end

