function [] = sequentialClassification(a,b)
    % 1. Learn three sequential classifers, and for each one plot the resulting
    % classifcation boundary along with the data points.
    for i = 1:3
        tic
        [sequentialDiscriminantClassifier, X, Y] = sequentialDiscriminantGenerator(a,b,999);
        toc
        figure;
        plotClasses(a,'Class A',b,'Class B');
        hold on;
        contour(X,Y,sequentialDiscriminantClassifier, 'DisplayName', 'SD Boundary');
        title(sprintf('Sequential Classifier %d',i));
    end

    % 3. In the above development we did not limit the number of sequential
    % classifers. Suppose we limit the sequential classifer to J classifers
    % G1; : : : ;GJ . We want to see how the experimental error rate varies
    % with J. For each value of J = 1; 2; : : : ; 5, learn a sequential classifer
    % 20 times to calculate the following:
    % (a) the average error rate
    % (b) minimum error rate
    % (c) maximum error rate
    % (d) standard deviation of the error rates
    % Produce a plot showing these results as a function of J.
    j_limit_range = 1:5;
    for j_limit = j_limit_range
        fprintf('J-Limit = %d\n', j_limit);
        tic
        for i = 1:20
            [sequentialDiscriminantClassifier, X, Y] = sequentialDiscriminantGenerator(a,b,j_limit);
            Sequential_classify = classifyPoints(X, Y, sequentialDiscriminantClassifier, a, 1, b, 2);
            conf_Sequential = confusionmat(Sequential_classify(:,1), Sequential_classify(:,2));
            error_rate(j_limit,i)=(conf_Sequential(1,2) + conf_Sequential(2,1))/(length([a;b]));
        end
        toc
        average_error_rate(j_limit) = mean(error_rate(j_limit,:));
        min_error_rate(j_limit) = min(error_rate(j_limit,:));
        max_error_rate(j_limit) = max(error_rate(j_limit,:));
        std_error_rate(j_limit) = std(error_rate(j_limit,:),1); %IS THIS THE CORRECT standard dev?
    end
    % Plot as function of J
    figure;
    hold on;
    errorbar(j_limit_range,average_error_rate,std_error_rate,'DisplayName','Avg error rate')
    % plot(j_limit_range, average_error_rate,'DisplayName','Avg error rate')
    plot(j_limit_range, min_error_rate,'DisplayName','Min error rate')
    plot(j_limit_range, max_error_rate,'DisplayName','Max error rate')
    plot(j_limit_range, std_error_rate,'DisplayName','Std dev error rate')
    legend;
end