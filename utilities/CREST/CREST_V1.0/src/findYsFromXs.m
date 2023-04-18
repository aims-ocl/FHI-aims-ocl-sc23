%% Extrapolate a polynomial on given x-y data from given x's
%% Take given x's and x-y data, return y-data after interpolation and
% averaging, for enhanced numerical stability.
% Input:
% Xs - the x positions to return y's for
% xyData - y vs x data
% averagingInterval - interval on which to extrapolate+average
function ys_at_Xs = findYsFromXs(Xs, xyData, averagingInterval)

Xs = reshape(Xs, numel(Xs), 1);
xData = xyData(:, 1);
yData = xyData(:, 2);

% Find potentials point-by-point; interpolate and average in order to
% smooth for greater stability
ys_at_Xs = zeros(size(Xs));
for i = 1 : numel(Xs)
    lowerXLim = Xs(i) - averagingInterval/2;
    upperXLim = Xs(i) + averagingInterval/2;
    % find indices around point
    intervalLeftIndex = find(xData <= lowerXLim, 1, 'last');
    intervalRightIndex = find(xData >= upperXLim, 1, 'first');
    % Throw error if given point is not included in given potential
    if isempty(intervalLeftIndex) || isempty(intervalRightIndex)
        err = MException('MATLAB:findYVals:PointsOutsideData', ...
            'Requested extrapolation from a point outside the range or too close to the edges of the data');
        throw(err)
    end
    
    if intervalRightIndex == intervalLeftIndex
        % Indices equal: limits are both equal to an x point. Take
        % potential at that point
        ys_at_Xs(i) = yData(intervalLeftIndex);
    elseif intervalRightIndex - intervalLeftIndex == 1
        % Indices 1 apart: averaging interval  << spacing of z data.
        %  In this case simply take the point along the line in between
        ys_at_Xs(i) = yData(intervalLeftIndex) + ...
            (yData(intervalRightIndex) - yData(intervalLeftIndex)) ./ (xData(intervalRightIndex) - xData(intervalLeftIndex)) .* ...
            (Xs(i) - xData(intervalLeftIndex));
    else
        intervalIndices = (intervalLeftIndex : intervalRightIndex)';
        
        % Interpolate for better stability vs. numerical noise
        yInterp = spline(xData(intervalIndices), yData(intervalIndices));
        
        % Evaluate mean of interpolation in interval
        ys_at_Xs(i) = integral(@(x)ppval(yInterp, x), lowerXLim, upperXLim) ./ averagingInterval;
    end
end

end


