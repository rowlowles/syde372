function ModelEstimation1D(data, actualDistribution, label, range)

[ans_var, ans_mu, ans_lambda, a_ans, b_ans] = parametricEstimator1D(data);

estimate_a_gauss = gaussian(ans_mu, ans_var^0.5);
estimate_a_exp = exponential(ans_lambda);
estimate_a_uniform = uniform(a_ans,b_ans);

PlotDistributions(estimate_a_gauss, actualDistribution,...
    range,...
    "Estimate of Dataset " + label + " using Gaussian Parametric Estimate",...
    'GaussianParam_' + label)
PlotDistributions(estimate_a_exp, actualDistribution,...
    range,...
     "Estimate of Dataset " + label + " using Exponential Parametric Estimate"...
     , 'ExpParam_' + label)
PlotDistributions(estimate_a_uniform, actualDistribution,...
    range,...
    "Estimate of Dataset " + label + " using Uniform Parametric Estimate",...
    'UniformParam_' + label)


end

