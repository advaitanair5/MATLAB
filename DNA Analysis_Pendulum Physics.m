%% Homework 3
%% Advaita Nair
%% UID: 206250044

%%Clear Cache
clear all
close all
clc
pause (0.5)

%%Input homework problem number
homework = input('Enter the example number (1-2):');
while homework ~= 1 && homework ~=2
        fprintf('Error: Unknown homework problem number.\n')
        homework = input('Enter the example number (1-2):');
end
switch(homework)

    %%Problem 1: Three Species Problem
    case 1

        %start timer
        tic;

        %define time stepping patterns
        tf = 12;
        dt = 0.0005;

        %calculate number of steps
        nt = ceil(tf/dt);
        
        %set initial conditions
        x0 = 2.00;
        y0 = 2.49;
        z0 = 1.50;
        t = 0;

        %define lotka-volterra coefficients
        alpha_x = 0.70;
        beta_x = 1.25;
        gamma_x = 0.45;

        alpha_y = 1.00;
        beta_y = 0.75;
        gamma_y = 1.25;

        alpha_z = 1.50;
        beta_z = 1.00;
        gamma_z = 1.00;

        %create results table
        results = zeros(nt, 4); 

        %start the for loop
        for i = 1:1:nt

            %Forward Euler Updating Method
            xp1 = x0 + dt*(alpha_x*x0*(1 - x0/20) - beta_x*x0*y0 - ...
                gamma_x*x0*z0);
            yp1 = y0 + dt*(alpha_y*y0*(1 - y0/25) - beta_y*y0*x0 - ...
                gamma_y*y0*z0);
            zp1 = z0 + dt*(alpha_z*z0*(1 - z0/30) - beta_z*z0*x0 - ...
                gamma_z*z0*y0);

            %store the values
            results(i, :) = [t, round(x0, 4), round(y0, 4), round(z0, 4)];

            %Update
            x0 = xp1;
            y0 = yp1;
            z0 = zp1;

            t = t + dt;
        end

        %end timer
        elapsed_time = toc;

        %print values
        resultsTable = array2table(results, 'VariableNames', { ...
            'Time', 'X', 'Y', 'Z'});
        disp(resultsTable)
        fprintf('\n The elapsed time is %f', elapsed_time);

    %Problem 2: The Pocket Change Problem
    case 2

        %create cash transaction ranges
        cash = linspace(0, 99, 100);

        %define value of coins
        q = 24; %quarter
        d = 12; %dime
        n = 4; %nickel
        p = 1; %penny

        %initialize coin values &total coins
        num_q = 0;
        num_d = 0;
        num_n = 0;
        num_p = 0;
        total_coins = 0;

        %set a loop for every cash transaction
        for i = 1:length(cash)
            amount = cash(i);
            %calculate most amount of quarters
            q_rmd = mod(amount, q);
            num_q = (amount - q_rmd)/q;
            %calculate most amount of dimes after max quarters
            d_rmd = mod(q_rmd, d);
            num_d = (q_rmd - d_rmd)/d;
            %calculate most amount of nickels after max dimes
            n_rmd = mod(d_rmd, n);
            num_n = (d_rmd - n_rmd)/n;
            %calculate most amount of pennies after max nickels
            num_p = n_rmd;
            %add up num of coins per transaction
            num_coins = num_q + num_d + num_n + num_p;
            %add up total num of coins for all transactions
            total_coins = total_coins + num_coins;
        end

        %calculate the average number of coins for all transactions
        avg_coin = total_coins/length(cash);

        %print
        fprintf('Average Number of Coins = %.2f', avg_coin);

end 

            

        