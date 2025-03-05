%Homework 5
%Advaita Nair
%UID: 206250044

function [ xa ] = averagePts( xs, w)

%initialize array
xa = zeros(1, length(xs));

%error statement
if sum(w) == 0
    error('Enter a non-zero value for the vector');
end

%weighted average
w = w/sum(w);

%create for loop
for i = 1:length(xs)
    %first element
    if i == 1
        im1 = length(xs);
    else
        im1 = i-1;
    end
    %last element
    if i == length(xs)
        ip1 = 1;
    else
        ip1 = i+1;
    end

    % Weighted Average of Positions
    xa(i) = w(1)*xs(im1) + w(2)*xs(i) + w(3)*xs(ip1);
end

