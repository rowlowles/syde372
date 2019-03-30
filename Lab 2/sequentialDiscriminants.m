function [] = sequentialDiscriminants(classA,classB)
    j = 1;
    [X, Y, classifier_MED] = MED(classA,classB);
    % calculate confusion matrix
    % Calculate error
    MED_classify = classifyPoints(X, Y, classifier_MED, classA, 1, classB, 2);
    conf_ged1 = confusionmat(MED_classify(:,1), MED_classify(:,2));
    
end

