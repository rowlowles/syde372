close all;
load('lab2_1');
load('lab2_2');

% Parametric Estination - 1D
[ans_sigma, ans_mu, ans_lambda] = parametricEstimator1D(a);

distribution_a_actual = gaussian(ans_mu, ans_lambda);


% Parametric Estimation - 2D
parametricEstimation2D(al,bl,cl);