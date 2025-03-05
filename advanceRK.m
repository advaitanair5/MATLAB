%Homework 5
%Advaita Nair
%UID: 206250044

%runge kutta function
function [ y ] = advanceRK(y, dt, method)
    
    % Define time stepping patterns
    t_half = 2.45;
    nt = length(y);
    decay_rate = -log(2) / t_half;

    % Runge-Kutta defined
    for k = 1:nt-1
        c1 = dt * decay_rate * y(k);
        c2 = dt * decay_rate * (y(k) + 0.5 * c1);
        c3 = dt * decay_rate * (y(k) + 0.5 * c2);
        c4 = dt * decay_rate * (y(k) + c3);
        
        % Check which order
        if method == 1
            y(k+1) = y(k) + c1;  % Forward Euler (RK1)
        elseif method == 2
            y(k+1) = y(k) + c2;  % RK2
        elseif method == 4
            y(k+1) = y(k) + (1/6)*c1 + (1/3)*c2 + (1/3)*c3 + (1/6)*c4;% RK4
        end
    end
end