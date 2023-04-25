close;
clear;
clc;

% Constants
P0  = 100;
sigma   = 0.15;
% mu  = 0.03;
mu  = 0.05*sigma^2/2;
delT    = 1/260;
days    = 260*50;

% Finding average return
num_sim     = 100000;
Pf_vec  = 1:num_sim;
for i = 1:num_sim
    [t_vec, P_vec] = stock_model(P0, mu, sigma, delT, days);
    Pf_vec(i) = P_vec(end);
end
disp("Average value after 20 years is: " + mean(Pf_vec));

% Plots
% for plt = 1:2
%     figure(plt);
%     clf;
%     hold on;
%     for i = 1:5
%         % Find new random stonks vector
%         [t_vec, P_vec] = stock_model(P0, mu, sigma, delT, days);
%         plot(t_vec, P_vec, 'LineWidth', 1);
%     end
%     xlabel("time (day)");
%     ylabel("Stock Price");
% end