function [w, mu, sigma, P0, names] = create_portfolio(filepath, alpha)
% Set the maximum proportion of a stock to buy
% MAX_PROP = 0.25;
MAX_PROP = 1;
% Set the number of trading days per year
DAYS_PER_YEAR = 252;
% Read our stock prices
stocks = readtable(filepath);
% Read our stock tickers
names = stocks.Properties.VariableNames(2:end);
% Load the raw stock history
P = table2array(stocks(:,2:end));
% Get the current prices
P0 = P(end,:);
% Get the number of stocks being tracked
[~, n] = size(P);
% Get the daily returns for each stock
R = (P(2:end, :)-P(1:end-1, :))./P(1:end-1, :);
% Get the average daily return for each stock
mu = mean(R);
% Get the covariance matrix between stocks
C = cov(R);
% Load the duration being used
times = table2array(stocks(:,1));
duration = between(times(1), times(end));
years = calyears(duration) + (calmonths(duration)/12) + (caldays(duration)/365);
delta_t = years * DAYS_PER_YEAR / n;
% Calculate our volatilities
sigma = zeros(1, n);
for i = 1 : n
    sigma(i) = sqrt(C(i,i))/delta_t;
end
% Solve the thing
[w, optVal] = quadprog((1-alpha)*2*C, -alpha*mu, [], [], ...
    ones(1, n), [1], zeros(n, 1), MAX_PROP*ones(n, 1));
% If it suggests investing less than 0.1% into a stock, call that zero
MIN_PERCENT = 0.001;
for i = 1 : n
    if w(i) < MIN_PERCENT
        w(i) = 0;
    end
end
% Re-scale our weights to account for any lost values
w = w / sum(w);
% Load our weights into a table
weights = cell2table(reshape(names, [n 1]), "VariableNames", ["Ticker"]);
weights.Weight = reshape(w, [n 1]);
% Sort by weight, then alphabetically
weights = sortrows(weights, 1);
weights = sortrows(weights, 2, "descend"); % REMOVE SEMICOLON TO VIEW TABLE
end
