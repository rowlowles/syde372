function plotStdContours(stdRange, varargin)
% varargin set of (mean, sigma) 
for index = 1:2:size(varargin,2)
    mean = varargin{index};
    sigma = varargin{index+1};
    for i = 1:size(stdRange,2)
        plotStdEllipse(sigma,mean,stdRange(i));
        hold on;
    end
end

end

