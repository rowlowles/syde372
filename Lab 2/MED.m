function [feature1Vals, feature2Vals, classifier_MED_1] = MED(classA,classB)
    muA = datasample(classA,1);
    muB = datasample(classB,1);

    minValue = floor(min(min(classA, classB)));
    maxValue = ceil(max(max(classA, classB)));

    feature1Vals = minValue(1):0.5:maxValue(1);
    feature2Vals = minValue(2):0.5:maxValue(2);

    [X, Y] = meshgrid(feature1Vals, feature2Vals);

    arrSize = [size(feature2Vals,2) size(feature1Vals,2)];

    classifier_MED_1 = zeros(arrSize);

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
    
    figure;
    plotClasses(classA,'Class A',classB,'Class B');
    hold on;
    contour(X,Y,classifier_MED_1, 'DisplayName', 'MED Boundary');
end

