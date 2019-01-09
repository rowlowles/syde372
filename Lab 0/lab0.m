clear all;
close all;

% Mean and variance of Probability Density Function
mu = [0; 0];
sigma = [1 0; 0 1];

stepSize = .1;

x1Range = [-3:stepSize:3];
x2Range = [-3:stepSize:3];

y = Gauss2d(x1Range, x2Range, mu, sigma);

% Show a 3-D plot of the pdf
figure
subplot(2,1,1);
surf(x1Range,x2Range,y);
title('Surface Plot')
xlabel('x_1');
ylabel('x_2');
% Show contours of the pdf
subplot(2,1,2);
contour(x1Range,x2Range,y);
title('Contour Plot')
xlabel('x_{1}');
ylabel('x_{2}');
axis equal
% Show a colour map of the pdf
figure
imagesc(x1Range,x2Range,y)
title('Value Plot')
xlabel('x_{1}');
ylabel('x_{2}');

z = (y>.1);

figure
imagesc(x1Range,x2Range,z)
hold on % allow us to plot more on the same figure
plot(mu(1,1),mu(2,1),'rx'); % plot the mean
xlabel('x_{1}');
ylabel('x_{2}');