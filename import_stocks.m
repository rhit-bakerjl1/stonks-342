function [history, names, growth_rate, volatility, covariance] = import_stocks(filepath)
% Parses the data in the inputted .csv file. Returns the stock price
% history, stock names, annual growth rates, and volatilities.

% Set days per year to 260 (approx. number of trading days)
DAYS_PER_YEAR = 260;
% Read stock prices from a .csv generated from scrape_stocks.py
stocks = readtable(filepath);
% Extract the stock price names
names = stocks.Properties.VariableNames(2:end);
% Get the duration of time being used
times = table2array(stocks(:,1));
duration = between(times(1), times(end));
years = calyears(duration) + (calmonths(duration)/12) + (caldays(duration)/365);
% Extract the stock price history into an array
history = table2array(stocks(:,2:end));
% Extract the size of our stock history
% Rows = number of days, cols = number of stocks
[rows, cols] = size(history);
% Get the time step delta_t
delta_t = years * DAYS_PER_YEAR / rows;
% Calculate the covariance of stocks from the price history
covariance = cov(history);
% Calculate growth rate and volatility of each stock
growth_rate = zeros(1, cols);
volatility = zeros(1, cols);
for i = 1 : cols
    growth_rate(i) = (history(end, i) - history(1, i))*years/history(1, i);
    volatility(i) = sqrt(covariance(i,i))/delta_t;
end