function parametricEstimation2D(a, b, c)

mean_a = mean(a);
cov_a = cov(a);

mean_b = mean(b);
cov_b = cov(b);

mean_c = mean(c);
cov_c = cov(c);

[feature1Vals, feature2Vals, MAPClassificationABC] = ML(a, mean_a, cov_a, b, mean_b, cov_b, c, mean_c, cov_c);

figure
hold on;
set(gca, 'ydir', 'normal');
%axis equal;
plotClasses(a,'Class A',b,'Class B', c,'Class C');
[X_ML, Y_ML] = meshgrid(feature1Vals, feature2Vals);
contour(X_ML,Y_ML,MAPClassificationABC,'DisplayName','ML boundary');


end

