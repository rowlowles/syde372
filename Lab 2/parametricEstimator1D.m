function parametricEstimator1D(input_data)
    syms mu sigma_x lambda a b x
    gauss = @(sigma_x, mu, x) (1/(2*pi*sigma_x^2)^(1/2)) * exp(-0.5*(x-mu)^2/(sigma_x^2));
    exponential = @(lambda, x) lambda*exp(-lambda*x);
%     uniform = @(a,b,x) piecewise((a<x)&&(x<b), 1/(b-a) * (b-a)/2, 0);
    uniform = @(a,b,x) (b+a)/2;

    
    log_likelihood_gauss = 0;
    log_likelihood_exponential = 0;
    likelihood_uniform = 1;
    for i = 1:length(input_data)
         log_likelihood_gauss = log_likelihood_gauss + log(gauss(sigma_x, mu, input_data(i)));
         log_likelihood_exponential = lqog_likelihood_exponential + log(exponential(lambda,input_data(i)));  
         likelihood_uniform = likelihood_uniform * uniform(a, b, input_data(i));
    end
%     dGauss_dSigma = diff(log_likelihood_gauss, sigma_x);
%     dGauss_dMu = diff(log_likelihood_gauss, mu);
%     ans_mu = double(solve(dGauss_dMu, mu));
%     ans_sigma = double(solve(subs(dGauss_dSigma,mu, ans_mu), sigma_x));
%     
%     dExp_dLambda = diff(log_likelihood_exponential, lambda);
%     ans_lambda = double(solve(dExp_dLambda, lambda));
    likelihood_uniform = log(likelihood_uniform);
    d_likelihood_uniform_d_a = diff(likelihood_uniform, a);
    d_likelihood_uniform_d_b = diff(likelihood_uniform, b);
    
    a_ans = (solve(d_likelihood_uniform_d_a, a));
    b_ans = solve(d_likelihood_uniform_d_b, b);
    
    a_ans = subs(a_ans(b, b_ans))

    
    std(input_data)
    mean(input_data)
end

