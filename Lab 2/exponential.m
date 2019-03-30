function [output] = exponential(lambda)
output = @(X)lambda*exp(-lambda*X);
end