function plotStdEllipse(sigma,mean, std, displayName)
    [V, D] = eig(sigma.*std);

    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];

    ellipsePlot = plot(a(1, :) + mean(1), a(2, :) + mean(2), 'DisplayName', displayName);
    ellipsePlot.LineWidth = 3;
    
end

