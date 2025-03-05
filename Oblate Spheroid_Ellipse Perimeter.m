%% Homework 1
%% Author: Advaita Nair 
%% UID: 206250044

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

        %%input r1 and r2
        r1 = input('Please enter a value for r1: \n');
        r2 = input('Please enter a value for r2: \n');
        
        %%error statement for which r1 & r2 are invalid
        if r1 >= r2
            error('r2 must be greater than r1.')
        elseif r1<=0||r2<=0
            error('The radii must be greater than 0.')
        end
        
        %Define formulas
        gamma = acos(r2/r1);
        area= 2*pi*(r1^2 +r2^2/sin(gamma)*log(cos(gamma)/(1-sin(gamma))));
        approx = 4*pi*((r1 + r2)/2)^2;
        
        %prints the values, rounded to 10 digits
        fprintf('\n');
        fprintf('For a radius of r1 = %f\n', r1);
        fprintf('and for a radius of r2 = %f\n', r2);
        fprintf('Area = %.10g\n', area);
        fprintf('Approx = %.10g\n', approx);

    %%Problem 2
    case 2

        %input a and b
        a = input('Please enter a value for a. \n');
        b = input('Please enter a value for b. \n');
        
        %error statement, a and/or b <=0 would not work physically
        if a<=0||b<=0
            error('a and/or b cannot be less than or equal to 0.')
        end
        
        %define formulas
        h = ((a - b)/(a + b))^2;
        p1 = pi*(a + b);
        p2 = pi*sqrt(2*(a^2 + b^2));
        p3 = pi*sqrt(2*(a^2 + b^2)-((a - b)^2)/2);
        p4 = pi*(a + b)*(1+(h/8))^2;
        p5 = pi*(a + b)*(1 + ((3*h)/(10 + sqrt(4 - 3*h))));
        p6 = pi*(a + b)*((64 - 3*h^2)/(64 - 16*h));
        p7 = pi*(a + b)*((256 - 48*h - 21*h^2)/(256 - 112*h + 3*h^2));
        p8 = pi*(a + b)*((3 - sqrt(1 - h))/2);
        
        %prints the values (perimeter (when a==b) defined in if statement)
        if a==b
            r=a;
            p = 2*pi*r;
            fprintf('\n');
            fprintf('For an entered value of a = %f\n', a);
            fprintf('For an entered value of b = %f\n', b);
            fprintf('h = %f\n', h);
            fprintf('P = %f\n', p);
            fprintf('P1 = %f\n', p1);
            fprintf('P2 = %f\n', p2);
            fprintf('P3 = %f\n', p3);
            fprintf('P4 = %f\n', p4);
            fprintf('P5 = %f\n', p5);
            fprintf('P6 = %f\n', p6);
            fprintf('P7 = %f\n', p7);
            fprintf('P8 = %f\n', p8);
        else
            fprintf('\n');
            fprintf('For an entered value of a = %f\n', a);
            fprintf('For an entered value of b = %f\n', b);
            fprintf('h = %f\n', h);
            fprintf('P1 = %f\n', p1);
            fprintf('P2 = %f\n', p2);
            fprintf('P3 = %f\n', p3);
            fprintf('P4 = %f\n', p4);
            fprintf('P5 = %f\n', p5);
            fprintf('P6 = %f\n', p6);
            fprintf('P7 = %f\n', p7);
            fprintf('P8 = %f\n', p8);
        end
end