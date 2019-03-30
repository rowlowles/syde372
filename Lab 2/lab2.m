close all;
load('lab2_1');
load('lab2_2');

%%
% Parametric Estination - 1D
a_distribution = gaussian(5,1);
b_distribution = exponential(1);

ModelEstimation1D(a,a_distribution, "A", 0:0.01:10);
ModelEstimation1D(b,b_distribution, "B", 0:0.01:10);

%%
% Parametric Estimation - 2D
parametricEstimation2D(al,bl,cl);

%%
load('lab2_3');
sequentialDiscriminants(a,b);