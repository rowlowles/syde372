function [class] = sequentialDiscriminantApplier(Master_classify, Master_n_a_B, Master_n_b_A, X, Y, point)
    j = 1;
    while(1)
        j_classification = findClassification(X,Y,Master_classify(:,:,j),point(1), point(2));
        if(j_classification == 2 && Master_n_a_B(j) == 0)
            class = 2;
            return;
        end
        if(j_classification == 1 && Master_n_b_A(j) == 0)
            class = 1;
            return;
        end
        j = j + 1;
    end
end