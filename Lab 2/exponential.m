function [output] = exponential(lambda)
    function y = outputFunction(X)
        y = lambda*exp(-lambda*X);
    end

output = @outputFunction; 
end


