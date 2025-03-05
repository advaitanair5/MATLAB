function [M1, x_real, W_gas] = Quasi_work(phi, Omega, sig, eps, alpha)
    M=9; 
    m=1.5; 
    r=0.05; 
    Patm=101325;

    T=310; 
    n=0.323; 
    R=8.314; 
    g=9.81; 
    l1=1; 
    l2=0.03;
    mu_irrev=1450;
    mu_rev=0;
    
    % VDM coefficient
    Tc = 33.19; %K
    Pc = 1.313*10^6; %pa
    omega = -0.216;
    Tr = T/Tc;

    % P = C1/(V-b) - C2/((V+eps*b)*(V+sig*b))
    a = phi*alpha*R^2*Tc^2/Pc; %a(T)
    b = Omega*R*Tc/Pc;
    C1 = R*T;
    C2 = a;
    
    %% Initial mass
    M1 = -m-(Patm-R*T/(pi*r^2*l1/n-b)+a/((pi*r^2*l1/n+eps*b)* ...
        (pi*r^2*l1/n+sig*b)))*pi*r^2/g;
    masses = M1:1e-2:9;

    %% Define equations
    syms x M
    eqn = (M - M1) + ...
        (1 / (pi * r^2 * l1 / n - b) - 1 / (pi * r^2 * x / n - b))*...
        pi * r^2 * R * T / g + ...
        (-a / ((pi * r^2 * l1 / n + eps * b) * ...
        (pi * r^2 * l1 / n + sig * b)) + a / ...
        ((pi * r^2 * x / n + eps * b)*(pi * r^2 * x/n + sig* b)))*...
        pi * r^2/g == 0;

    P = C1/(pi*r^2*x/n-b)-C2/((pi*r^2*x/n+eps*b)*(pi*r^2*x/n+sig*b));

    x_real = zeros(size(masses));  % Preallocate

    for i = 1:length(masses)
        M_val = masses(i);
        eqn_M = subs(eqn, M, M_val);
        x_sol = vpasolve(eqn_M, x, l1);
        x_filtered = double(x_sol(imag(x_sol) == 0 & x_sol > 0.1));  % Filter real and positive solutions

        if ~isempty(x_filtered)  % Only store valid solutions
            x_real(i) = x_filtered(1);  % Take the first valid solution
            valid_indices(i) = true;
        end
    end

    % Filter out invalid values
    x_real = double(x_real(valid_indices));
    fprintf('Initial mass = %.6f\n', M1);

    %% Work done by gas
    W = int(P*pi*r^2,x,l1,x);
    W_gas_function = matlabFunction(W);
    W_gas = W_gas_function(x_real);
end