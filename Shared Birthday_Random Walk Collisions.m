%%Homework 6
%%Advaita Nair
%%UID: 206250044

%The Shared Birthday Problem: The script is designed to create random bdays
%and compare them to see if any birthdays are within the same week
%
%Random Walk Collisions: The script is designed to model random movement of
%particles and predict when they will collide


%Homework 6 Code:

%%Clear Cache
clear all
close all
clc
pause (0.5)

%%Question Switch

%%Input homework problem number
homework = input('Enter the example number (1-2):');
while homework ~= 1 && homework ~=2
        fprintf('Error: Unknown homework problem number.\n')
        homework = input('Enter the example number (1-2):');
end
switch(homework)

    %%Problem 1
    case 1
        %initialize variables
        num_trials = 1e4;
        num_req = zeros(1, num_trials);
        days_year = 365;
        rng('default');

        %a for loop that will run for the set amount of trials
        for trial = 1:num_trials
            %set initial conditions for each trial
            birthdays = [];
            match = false;
            n = 0;

            %create a while loop until a match is found
            while ~match
                %update variables 
                n = n+1;
                new_birthday = randi(days_year);
                birthdays = [birthdays, new_birthday];
                
                    %check matches for all previous birthays in array
                    for i=1:length(birthdays)-1
                        %calculate the diff for each prev bday to new one
                        day_diff = abs(birthdays(i)-new_birthday);

                        %check if it is within the same week & include
                        %wrapping around
                        if (day_diff < 7) || (days_year - day_diff < 7)
                            match = true;
                            break; %exits loop when match is found
                        end
                    end
            end
            %record the number of ppl until a match is found
            num_req(trial) = n;
        end
        
        %plot histogram of n values
        histogram(num_req);
        ylabel('Frequency')
        xlabel('Number of People Until Match')

        %calculate the median number of ppl before a match is found
        med_ppl = median(num_req);
        fprintf('Median Number of People = %02d\n', med_ppl)

                  

    %%Problem 2
    case 2
        %set seed
        rng('default');

        %set particle A initial position
        xA = -5;    xAk = xA;
        yA = 0;     yAk = yA;

        %set particle B initial position
        xB = 5;     xBk = xB;
        yB = 0;     yBk = yB;

        %set boundary conditions
        BC = [5, -5, -5, 5];

        %initialize variables
        n = 1000;    n0_vals = [];
        max_trial = 5000;   trial = 0;

        %create while loop for set max trials
        while trial <= max_trial

            %reset variables
            n0 = 0;
            collision = 0;
            xAk = -5;
            yAk = 0;
            xBk = 5;
            yBk = 0;

            %for one run
            while collision == 0 && n0<n
    
                %decide next position
                [xAkp1, yAkp1] = RandWalk(xAk, yAk, BC);
                [xBkp1, yBkp1] = RandWalk(xBk, yBk, BC);
    
                % Create Particle A on Grid for Step (k)
                xAk_val = [xAk - 0.5, xAk + 0.5, xAk + 0.5, xAk - 0.5];
                yAk_val = [yAk - 0.5, yAk - 0.5, yAk + 0.5, yAk + 0.5];
                % Create Particle A on Grid for Step (k+1)
                xAkp1_val = [xAkp1 - 0.5, xAkp1+0.5, xAkp1+0.5, xAkp1-0.5];
                yAkp1_val = [yAkp1 - 0.5, yAkp1-0.5, yAkp1+0.5, yAkp1+0.5];
                % Create Particle B on Grid for Step (k)
                xBk_val = [xBk - 0.5, xBk + 0.5, xBk + 0.5, xBk - 0.5];
                yBk_val = [yBk - 0.5, yBk - 0.5, yBk + 0.5, yBk + 0.5];
                % Create Particle B on Grid for Step (k+1)
                xBkp1_val = [xBkp1 - 0.5, xBkp1+0.5, xBkp1+0.5, xBkp1-0.5];
                yBkp1_val = [yBkp1 - 0.5, yBkp1-0.5, yBkp1+0.5, yBkp1+0.5];
    
                %plot figure for the first trial
                if trial == 0 
                    figure(1)
                    hold on
                    set(gca,'xtick',-5.5:1:5.5)
                    set(gca,'ytick',-5.5:1:5.5)
                    grid on
                    xlim([-5.5 5.5])
                    ylim([-5.5 5.5])
                    axis square
                    fill(xAk_val,yAk_val,'r')
                    fill(xAkp1_val,yAkp1_val,'b')
                    fill(xBk_val,yBk_val,'y')
                    fill(xBkp1_val,yBkp1_val,'g')
                    title('2D Random Walk','FontSize',24)
                    set(gcf,'Position',[30 350 600 600])
                    set(gca,'LineWidth',3,'FontSize',20)
                    hold off
                end

                %update new position
                xAk = xAkp1;    yAk = yAkp1;
                xBk = xBkp1;    yBk = yBkp1;
                n0 = n0 + 1;
    
                %determine if collision happened
                if xAk == xBk && yAk == yBk
                    collision = 1;
                end
            end

            %store num of steps before collision
            n0_vals = [n0_vals, n0];
            %update trials
            trial = trial + 1;
        end

        %print median
        med_n = median(n0_vals);
        fprintf('Median %.2f\n', med_n);

end