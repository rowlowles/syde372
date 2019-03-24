function parametricEstimator1D(input_data)
    syms mu sigma_x lambda a b x
    gauss = @(sigma_x, mu, x) 1/(2*pi*sigma_x^2)^1/2 * exp(-0.5*(x-mu)^2/(sigma_x^2));
    exponential = @(lambda, x) lambda*exp(-lambda*x);
    uniform = @(a,b,x) piecewise((a<x)&&(x<b), 1/(b-a) * (b-a)/2, 0);
    
    log_likelihood_gauss = 0;
    for i = 1:length(input_data)
       log_likelihood_gauss = log_likelihood_gauss + log(gauss(sigma_x, mu, input_data(i)));
%        log_likelihood_exponential = log_likelihood_exponential + exponential(lambda, 
        
    end
    dGauss_dSigma = diff(log_likelihood_gauss, sigma_x)
    dGauss_dMu = diff(log_likelihood_gauss, mu)
    ans_mu = double(solve(dGauss_dMu, mu))
%     ans_sigma = solve(subs(dGauss_dSigma, ans_mu), sigma_x)

    cov(input_data)
    mean(input_data)
end

