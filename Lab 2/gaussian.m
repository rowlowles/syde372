function [output] = gaussian(mu, sigma)
    output = @(X)(1/(2*pi*sigma)^(1/2)) * exp(-0.5*(X-mu).^2/(sigma)); 
end

