function plotStdContours(stdRange, varargin)
% varargin set of (mean, sigma, name) 
legend('-DynamicLegend');
for index = 1:3:size(varargin,2)
    mean = varargin{index};
    sigma = varargin{index+1};
    name = varargin{index+2};
    for i = 1:size(stdRange,2)
        plotStdEllipse(sigma,mean,stdRange(i),name);
        hold on
    end
end

end

