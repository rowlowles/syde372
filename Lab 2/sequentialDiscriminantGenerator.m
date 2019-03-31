function [Master_Master_classify, X, Y] = sequentialDiscriminantGenerator(classA,classB,j_limit)
    j = 1;
    minValue = floor(min(min(classA, classB)));
    maxValue = ceil(max(max(classA, classB)));
    X = minValue(1):1:maxValue(1);
    Y = minValue(2):1:maxValue(2);
    while (j<=j_limit)
        [classifier_MED] = MED(classA,classB, X,Y);
        % calculate confusion matrix
        % Calculate error
        MED_classify = classifyPoints(X, Y, classifier_MED, classA, 1, classB, 2);
        conf_MED = confusionmat(MED_classify(:,1), MED_classify(:,2));
        n_a_B = conf_MED(1,2);
        n_b_A = conf_MED(2,1);
        if(n_a_B ~= 0 && n_b_A ~= 0)
            continue %try again
        end
        Master_classify(:,:,j) = classifier_MED;
        Master_n_a_B(j) = n_a_B;
        Master_n_b_A(j) = n_b_A;
        
        if(n_a_B == 0 && ~isempty(classB))
            % remove these points
            indexes_b_B = find(MED_classify(:,1) == 2 & MED_classify(:,2) == 2) - size(classA,1);
            classB(indexes_b_B,:) = [];
        end
        if(n_b_A == 0 && ~isempty(classA))
            indexes_a_A = find(MED_classify(:,1) == 1 & MED_classify(:,2) == 1);
            classA(indexes_a_A,:) = [];
        end
        if(isempty(classA) && isempty(classB))
            break;
        end
        
        j = j + 1;
    end
    
%     Master_Master_classify = zeros(size(Master_classify(:,:,1)));
%     for i = 1:length(X)
%         for j = 1:length(Y)
%             Master_Master_classify(j,i) = sequentialDiscriminantApplier(Master_classify, Master_n_a_B, Master_n_b_A, X, Y, [X(i),Y(j)]);
%         end
%     end
    Master_Master_classify = zeros(size(Master_classify(:,:,1)));
    for i = 1:size(Master_classify,3)
        if(Master_n_a_B(i) == 0)
            %class = 2
            Master_Master_classify((Master_Master_classify==0)&(Master_classify(:,:,i)==2)) = 2;
        end
        if(Master_n_b_A(i) == 0)
            %class = 1
            Master_Master_classify((Master_Master_classify==0)&(Master_classify(:,:,i)==1)) = 1;
        end
        if(i==size(Master_classify,3))
            Master_Master_classify((Master_Master_classify==0)&(Master_classify(:,:,i)==2)) = 2;
            Master_Master_classify((Master_Master_classify==0)&(Master_classify(:,:,i)==1)) = 1;
        end
    end
end

