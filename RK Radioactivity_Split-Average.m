%% Homework 4
% Advaita Nair
% UID: 206250044
%% The Split-and-Average Problem: 
% The script is designed to find the midpoint of every-ther element and 
% find the weighted average of the 3 nearest neighbors.
%% Runge-Kutta Radioactivity:
% The script is designed to find the amount of carbon remaining using 
% 1st, 2nd, and 4th order Runge-Kutta and compare that to the exact formula
% using different time steps for comparison.

%% The Script

%%Clear Cache
clear all
close all
clc
pause(0.5)

%%Input homework problem number
homework = input('Enter the example number (1-2):');
while homework ~= 1 && homework ~=2
        fprintf('Error: Unknown homework problem number.\n')
        homework = input('Enter the example number (1-2):');
end
switch(homework)

    %%Problem 1: The Split-and-Average Problem
    case 1

        %array to be tested
        x = [0, 0, 1, 1];
        y = [0, 1, 0, 1];

        figure;
        plot(x, y, 'bo');
        hold on;

        %set conditions 
        maxdisp = 1e-3;
        iteration = 0;
        maxiteration = 100;

        maxdxdy = Inf;

        while maxdxdy > maxdisp && iteration < maxiteration
           
            %call split function
            xs = splitPts(x);
            ys = splitPts(y);

            %call average function
            w = [.25, .25, .25];
            xa = averagePts(xs, w);
            ya = averagePts(ys, w);

            %compute displacement between new & old pt
            dx = xa - xs;
            dy = ya - ys;
            disp = sqrt(dx.^2 + dy.^2);
            maxdxdy = max(disp);

            %update for next iteration
            x = xa;
            y = ya;
            iteration = iteration + 1;
        end

        plot(x, y, 'r*');
        title('Initial and Final Point Distribution');
        legend('Initial', 'Final');
        hold off;


    %%Problem 2: Runge-Kutta Radioactivity
    case 2

        %rk orders
        methods = [1, 2, 4];

        %initialize
        t0 = 0;
        tf = 15;
        dt_vals = [1, 0.1, 0.01];
        t_half = 2.45;
        errors = zeros(length(dt_vals), length(methods));
     
        %call function for each diff dt & method
        for i = 1:length(dt_vals)
            dt = dt_vals(i);

            nt = ceil(tf/dt);
            t = t0:dt:tf;

            %exact formula
            y = exp((-log(2)/t_half)*t);

            %create figure & hold on so multiple line can be plotted
            figure;
            hold on;

            %initialize legend array
            legends = cell(1, 4);

            %loop for each RK method
            for j = 1:length(methods)
                method = methods(j);

                %call RK function
                y_rk = zeros(1,nt); y_rk(1) = 1;
                y_rk = advanceRK(y, dt, method);

                %plot results
                plot(t, y_rk);
                
                %update legends array
                legends{j} = sprintf('RK Order %d', method);

                %compute average error
                avg_error = mean(abs(y_rk - y));
                errors(i, j) = avg_error;
           
            end

            %plot exact solution
            plot(t, y);
            legends{end} = 'Exact Solution';
            xlabel('Time Elapsed (s)');
            ylabel('Amount of Carbon Remaining');
            title(sprintf('Half Life of Carbon for dt = %.3f', dt));
            
            %display legend
            legend(legends);
           
            %stop holding onto lines
            hold off;
            grid on;
        end

        %print table of errors
        fprintf('\n%-10s %-10s %-10s %-10s\n', 'dt', ...
            'RK1', 'RK2', 'RK4');
        for i = 1:length(dt_vals)
            fprintf('%-10.3f: %-10.2e %-10.2e %-10.2e\n', dt_vals(i), ...
                errors(i, 1), errors(i, 2), errors(i, 3));
        end
end