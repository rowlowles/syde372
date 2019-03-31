function [classifier_MED_1] = MED(classA,classB,feature1Vals, feature2Vals)
    muA = datasample(classA,1,1);
    muB = datasample(classB,1,1);

    [X, Y] = meshgrid(feature1Vals, feature2Vals);

    arrSize = [size(feature2Vals,2) size(feature1Vals,2)];

    classifier_MED_1 = zeros(arrSize);
    boundary1 = [];

    for x1MED = 1:arrSize(2)
        for y1MED = 1:arrSize(1)
            pointCord = [feature1Vals(x1MED), feature2Vals(y1MED)];
            distanceA = sum((muA-pointCord).^2)^0.5;
            distanceB = sum((muB-pointCord).^2)^0.5;
            if distanceA < distanceB
                classifier_MED_1(y1MED,x1MED) = 1;
            else
                classifier_MED_1(y1MED,x1MED) = 2;
            end
        end
    end
    
%     figure;
%     plotClasses(classA,'Class A',classB,'Class B');
%     hold on;
%     contour(X,Y,classifier_MED_1, 'DisplayName', 'MED Boundary');

%     errorClass = size(find(attachedMat(:,1) ~= attachedMat(:,2)),1)/size(attachedMat,1);
end

