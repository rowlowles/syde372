function plotStdEllipse(sigma,mean, std)
    [V, D] = eig(normalize(sigma).*std);

    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];

    plot(a(1, :) + mean(1), a(2, :) + mean(2));
end

