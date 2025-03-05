%% Homework 4
%% Advaita Nair
%% UID: 206250044

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

    %%Problem 1: DNA Analysis
    case 1

        %load in the data file and convert dna to string 
        load('chr1_sect.mat');
        dna = num2str(dna);
        %dna = '143122413143122221413';
        numBases = length(dna);

        %initialize variables
        start = 0;
        stop = 0;
        protein_length = 0;
        protein_lengths = [];

        %create for loop that will check every codon (3 bases)
        for k = 1:3:numBases
            %Starting Codon
            if start == 0
                if dna(k)=='1' && dna(k+1)=='4' && dna(k+2)=='3'
                    start = k;
                end
            elseif start ~=0
            %Stop Codons
                if (dna(k)=='4' && dna(k+1)=='1' && dna(k+2)=='1')||...
                   (dna(k)=='4' && dna(k+1)=='1' && dna(k+2)=='3')||...
                   (dna(k)=='4' && dna(k+1)=='3' && dna(k+2)=='1')
                   stop = k+2;
                   %Calculate the number of proteins & append to array
                   protein_length = (stop - start+1);
                   protein_lengths = [protein_lengths, protein_length];
                   %reset start and stop codons for next protein
                   start = 0;
                   stop = 0;
                end
            end
        end
        %Analyze Results
        total_protein = numel(protein_lengths);
        avg_length = mean(protein_lengths);
        med_length = median(protein_lengths);
        max_length = max(protein_lengths);
        min_length = min(protein_lengths);

        %Print Out Results
        fprintf('Total Protein-Coding Segments: %.0f\n', total_protein);
        fprintf('Average Length: %.2f\n', avg_length);
        fprintf('Median Length: %.0f\n', med_length);
        fprintf('Maximum Length: %.0f\n', max_length);
        fprintf('Minimum Length: %.0f\n', min_length);

    %%Problem 2: The Pendulum Physics Problem
    case 2
        %Initial Set-Up
        g = 9.81;
        L = 2.5;
        theta_0 = pi/2;

        %define time stepping patterns
        t0 = 0;
        tf = 20;
        dt = 0.05;

        %calculate number of steps
        nt = ceil(tf/dt);
        t = linspace(t0,tf,nt);

        %preallocate array & assign initial values
        theta_e = zeros(1,nt); theta_e(1) = theta_0;
        omega_e = zeros(1,nt); omega_e(1) = 0; %starting from rest
        alpha_e = zeros(1,nt); alpha_e(1) = -g/L; %sin = pi/2 = 1
        energy_e = zeros(1,nt); energy_e(1) = g*L;

        %Forward Euler (EXPLICIT) Updating Method (for position, velocity, 
        % and acceleration)
        for k = 1:1:nt-1
            omega_e(k+1) = omega_e(k) + dt*(-g/L*sin(theta_e(k)));
            theta_e(k+1) = theta_e(k) + dt*omega_e(k);
            alpha_e(k+1) = -g/L*sin(theta_e(k));
            energy_e(k+1) = g*(L-L*cos(theta_e(k))) + 0.5*(L*omega_e(k))^2;
        end

        %plot all graphs (disp, vel, accel) & etot
        figure;
        hold on;
        plot(t, theta_e)
        plot(t, omega_e);
        plot(t, alpha_e);
        hold off;
        xlabel('Time Elapsed (s)');
        ylabel('Displacement(m), Velocity(m/s), Acceleration(m/s^2)')
        title('Displacement, Velocity, Acceleration vs. Time (EXPLICIT)');
        legend({'Displacement', 'Velocity', 'Acceleration'})
        grid on;
        
        figure;
        plot(t, energy_e);
        ylabel('Energy (J)');
        xlabel('Time Elapsed (s)');
        title('Total Energy vs Time (EXPLICIT)');
        grid on;
        
        %repeat for implicit
        theta_i = zeros(1,nt); theta_i(1) = theta_0;
        omega_i = zeros(1,nt); omega_i(1) = 0;
        alpha_i = zeros(1,nt); alpha_i(1) = -g/L*sin(theta_i(1));
        energy_i = zeros(1,nt); energy_i(1) = g*L;

        %Forward Euler Updating (IMPLICIT) Method (for position, velocity, 
        % and acceleration)
        for k = 1:1:nt-1
            omega_i(k+1) = omega_i(k) + dt*(-g/L*sin(theta_i(k)));
            theta_i(k+1) = theta_i(k) + dt*omega_i(k+1);
            alpha_i(k+1) = -g/L*sin(theta_i(k));
            energy_i(k+1)=g*(L-L*cos(theta_i(k+1)))+0.5*(L*omega_i(k+1))^2;
        end

        %plot all graphs (disp, vel, accel) & etot
        figure;
        hold on;
        plot(t, theta_i);
        plot(t, omega_i);
        plot(t, alpha_i);
        hold off;
        xlabel('Time Elapsed (s)');
        ylabel('Displacement(m), Velocity(m/s), Acceleration(m/s^2)')
        title('Displacement, Velocity, Acceleration vs. Time (IMPLICIT)');
        legend({'Displacement','Velocity','Acceleration'})
        grid on;
        
        figure;
        plot(t, energy_i);
        ylabel('Energy (J)');
        xlabel('Time Elapsed (s)');
        title('Total Energy vs Time (IMPLICIT)');
        grid on;
end