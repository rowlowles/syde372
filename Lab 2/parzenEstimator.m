function parzenEstimator( input_data, sigma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dims = size(input_data);
length = dims(2);
pdf = zeros(1,length);
h = 1/sqrt(length);
for j = 0.1:.1:10
    foo = int8(j*10);
    x = j;
    pointSum = 0;
    for i = 1:length
        
        u = -(x - input_data(i))^2/(2*sigma^2);
        gaussKernel = (1/(sigma*sqrt(2 * pi)))*exp(u);
        
        pointSum = pointSum + gaussKernel;
    end
    % j is x step, pointSum/length is y value
    pdf(foo) = pointSum/length;
end
% plot(pdf)
% x = pdf;
figure;
hold on;
titleString = sprintf('Parzen Estimation of Dataset A with Sigma = %.1d', sigma);
title(titleString)
plot(input_data)
plot(pdf)
legend('Actual','Parzen')

end

