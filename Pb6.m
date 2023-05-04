close;
clear;
clc;

% Setup for create_portfolio
N_alpha     = 20;
N_stocks    = 41;
alphas      = linspace(0,1,N_alpha);
filepath    = "stock_scraper/combined_history.csv";
days        = 21;
delT        = 1/260;

% Create Portfolio
weights_vec     = zeros(N_stocks, N_alpha);
[weights(:,1), mus, sigmas, P0_vec]  = create_portfolio(filepath, alphas(1));
for i = 2:length(N_alpha)
    [weights(:,i), ~, ~] = create_portfolio(filepath, alphas(i));
end

% Analyze my weights
V       = 2000;     % Dollars we invest
P0_mat  = P0_vec*ones(1,N_alpha);
x_mat   = weights*V./P0_mat;

