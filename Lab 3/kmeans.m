function [classifier] = kmeans(K, data)
    prototypes = datasample(data, K);
    prev_prototypes = zeros(size(prototypes));
    iter = 1;
    
    while ((iter < 10000) && ~(isequal(prev_prototypes,prototypes)))
        prev_prototypes = prototypes;
        classifier = zeros(length(data),1);
        for point = 1:length(data)
           closest=1;
           min_dist = inf;
           for proto = 1:K
               dist = sum((prototypes(proto,:)-(data(point,:))).^2);
               if(dist < min_dist)
                   min_dist = dist;
                   closest = proto;
               end
           end
           classifier(point) = closest;
        end

        for proto = 1:K
           prototypes(proto,:) = mean(data(classifier == proto,:)); 
        end
        iter = iter + 1;
    end
    
    kmean_fig = figure('units','normalized','outerposition',[0 0 0.75 0.75]);
    legend_data = cell(K*2,1);
    for classIdx = 1:K
        randcolor = [rand rand rand];
        current_data = data(classifier == classIdx,:);
        scatter(current_data(:,1),current_data(:,2),65,randcolor,'filled');
        hold on; 
        scatter(prototypes(classIdx,1),prototypes(classIdx,2),200,...
            randcolor,'*');
        hold on; 
        legend_data{classIdx*2-1} = sprintf('Class %d',classIdx);
        legend_data{classIdx*2} = sprintf('Mean for class %d',classIdx);
    end
    
    xlabel('Feature 1');
    ylabel('Feature 2');
    title('K-Means Classification for 32x32 Image Data');
    legend(legend_data, 'Location', 'northwest');
    saveas(kmean_fig,"Kmean_plot.png");
end

