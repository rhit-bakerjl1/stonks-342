function [sse] = sse(vec)
% Returns the sum of squared error from a simple bias regressor for
% an estimator.
b = ones(size(vec)) * mean(vec);
sse = sum((vec-b).^2);
end

