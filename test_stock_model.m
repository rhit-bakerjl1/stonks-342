close;
clear;
clc;

% Constants
P0  = 100;
mu  = 0.03;
sigma   = 0.15;
delT    = 1/260;
days    = 260*3;

% Plotting and finding new stuff
for plt = 1:2
    figure(plt);
    clf;
    hold on;
    for i = 1:5
        % Find new random stonks vector
        [t_vec, P_vec] = stock_model(P0, mu, sigma, delT, days);
        plot(t_vec, P_vec, 'LineWidth', 1);
    end
    xlabel("time (day)");
    ylabel("Stock Price");
end