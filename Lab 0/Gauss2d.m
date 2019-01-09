%GAUSS2D 2-D Gaussian pdf
%
% INPUTS:
% x1 - 1xN vector - range of x1 component
% x2 - 1xM vector - range of x2 component
% Mu - 2x1 vector - the mean of the random variable
% Sigma - 2x2 matrix - the variance matrix of
% the random variable
%
% OUTPUT:
% Y - NxM matrix - the pdf for the range of x1, x2
%
% USAGE (Example):
% Y = Gauss2d([-3:0.1:3],[0:0.1:6],[0 3]',[1 1; 1 4])
% This function creates the pdf for a 2D Gaussian
% random variable
function Y = Gauss2d(x1, x2, Mu, Sigma)

Y = zeros(length(x1), length(x2));
for i = 1:length(x1)
    for j = 1:length(x2)
        v = [x1(i); x2(j)];
        Y(i,j) = 1/(sqrt(2*pi*det(Sigma)))*exp(-.5*(v-Mu)'*inv(Sigma)*(v-Mu));
    end
end
        