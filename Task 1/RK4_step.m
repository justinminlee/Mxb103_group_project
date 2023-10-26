%% Function for RK4 Step
function [y_next, v_next] = RK4_step(y, v, C, K, L, g, dt)
    % K1
    K1_y = v;
    K1_v = g - C*abs(v)*v - max(0, K*(y - L));
    
    % K2
    K2_y = v + 0.5*dt*K1_v;
    K2_v = g - C*abs(K2_y)*K2_y - max(0, K*(y + 0.5*dt*K1_y - L));
    
    % K3
    K3_y = v + 0.5*dt*K2_v;
    K3_v = g - C*abs(K3_y)*K3_y - max(0, K*(y + 0.5*dt*K2_y - L));
    
    % K4
    K4_y = v + dt*K3_v;
    K4_v = g - C*abs(K4_y)*K4_y - max(0, K*(y + dt*K3_y - L));
    
    % Update
    y_next = y + (dt/6)*(K1_y + 2*K2_y + 2*K3_y + K4_y);
    v_next = v + (dt/6)*(K1_v + 2*K2_v + 2*K3_v + K4_v);
end
