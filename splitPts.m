%Homework 5
%Advaita Nair
%UID: 206250044


function [ xs ] = splitPts( x )
    %double array size to include split points
    xs = zeros(1, 2*length(x));
    xs(1:2:end) = x;

    %create a for loop to split the values
    for i = 1:1:length(x)
        if i ==length(x) %wrap around to first point
            xs(2*i) = (x(i) + x(1))/2;
        else
            xs(2*i) = (x(i) + x(i+1))/2;
        end
    end
end
