function [output] = uniform(a,b)
output = @(x) ((a<x).*(x<b))*(1/(b-a));
end

