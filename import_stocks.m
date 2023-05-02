function [history, names, growth_rate, volatility, covariance] = import_stocks(filepath)
% Parses the data in the inputted .csv file. Returns the stock price
% history, stock names, annual growth rates, and volatilities.

% Read stock prices from a .csv generated from scrape_stocks.py
stocks = readtable(filepath);
% Extract the stock price names
names = stocks.Properties.VariableNames(2:end);
times = table2array(stocks(:,1));
duration = between(times(1), times(end));
years = calyears(duration) + (calmonths(duration)/12) + (caldays(duration)/365);
history = table2array(stocks(:,2:end));
covariance = cov(history);
[~, cols] = size(history);
growth_rate = zeros(1, cols);
volatility = zeros(1, cols);
for i = 1 : cols
    growth_rate(i) = (history(end, i) - history(1, i))*years/history(1, i);
    volatility(i) = sqrt(sse(history(:,i)));
end