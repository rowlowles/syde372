function plotClasses(varargin)
% Plots classes 
    for i = 1:2:size(varargin,2)
        class = varargin{i};
        scatter(class(:,1),class(:,2));
        hold on; 
    end
    
    xlabel('Feature 1');
    ylabel('Feature 2');
    legend(varargin{2:2:size(varargin,2)})
end

