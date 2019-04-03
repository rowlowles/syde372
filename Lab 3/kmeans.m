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
    
    figure;
    for classIdx = 1:K
        current_data = data(classifier == classIdx,:);
        scatter(current_data(:,1),current_data(:,2),50,'filled');
        hold on; 
    end
    
    xlabel('Feature 1');
    ylabel('Feature 2');
    legend();
end

