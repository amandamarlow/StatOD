function [outputArg1,outputArg2] = householder(R_ap, b_ap, H, y)
%HOUSEHOLDER Summary of this function goes here
%   Detailed explanation goes here

A = [R_ap, b_ap; H, y];
[m, n] = size(H);
for k = 1:n
    
    for j = k+1:n+1
        
    end
end
end

