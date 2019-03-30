function [ans_var, ans_muVar, ans_lambda, a_ans, b_ans] = parametricEstimator1D(input_data)
    syms muVar var_x lambda a b x
    gauss = @(var_x, muVar, x) (1/(2*pi*var_x)^(1/2)) * exp(-0.5*(x-muVar)^2/(var_x));
    exponential = @(lambda, x) lambda*exp(-lambda*x);
    uniform = @(a,b,x) 1/(b-a); %assumes all training data comes from dist

    log_likelihood_gauss = 0;
    log_likelihood_exponential = 0;
    log_likelihood_uniform = 0;
    for i = 1:length(input_data)
         log_likelihood_gauss = log_likelihood_gauss + log(gauss(var_x, muVar, input_data(i)));
         log_likelihood_exponential = log_likelihood_exponential + log(exponential(lambda,input_data(i)));  
         log_likelihood_uniform = log_likelihood_uniform + log(uniform(a, b, input_data(i)));
    end
    
    dGauss_dSigma = diff(log_likelihood_gauss, var_x);
    dGauss_dmuVar = diff(log_likelihood_gauss, muVar);
    ans_muVar = double(solve(dGauss_dmuVar, muVar));
    ans_var = double(solve(subs(dGauss_dSigma,muVar, ans_muVar), var_x));
     
    dExp_dLambda = diff(log_likelihood_exponential, lambda);
    ans_lambda = double(solve(dExp_dLambda, lambda));
    
    d_likelihood_uniform_d_a = diff(log_likelihood_uniform, a);
    d_likelihood_uniform_d_b = diff(log_likelihood_uniform, b);
    %Note: these derivatives never reach zero. Also,
    %log_likelihood_uniform(a,b,x) is unbounded on 0<b-a<inf. It is largest
    %closest to b-a=0. Thus, choose b and a such that their difference is
    %as small as possible. To do so, choose the min training value for a,
    %and the largest training value for b. This satisfies the constraint
    %that all training data is within the estimated distribution, and that
    %b > a.
    a_ans = min(input_data);
    b_ans = max(input_data);

end