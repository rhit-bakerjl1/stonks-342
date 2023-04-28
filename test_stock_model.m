close;
clear;
clc;

% Constants
P0  = 100;
% sigma   = 0.15;
sigma = 1;
% mu  = 0.03;
mu  = 0.2*sigma^2/2;
% mu  = 0;
delT    = 1/260;
days    = 260*50;

% Finding average return
num_sim     = 10000;
Pf_vec  = 1:num_sim;
for i = 1:num_sim
    [t_vec, P_vec] = stock_model(P0, mu, sigma, delT, days);
    Pf_vec(i) = P_vec(end);
end
disp("Min value after 20 years is: " + min(Pf_vec));
disp("Max value after 20 years is: " + max(Pf_vec));
disp("Average value after 20 years is: " + mean(Pf_vec));
disp("Standard deviation after 20 years is: " + std(Pf_vec));

% Plots
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