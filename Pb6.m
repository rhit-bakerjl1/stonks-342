%% Analyzing Portfolio Choices
close;
clear;
clc;

% Setup for create_portfolio
N_alpha     = 100;
N_stocks    = 41;
N_iter      = 100000;
alpha_min   = 0;
alpha_max   = 1;
alphas      = linspace(alpha_min,alpha_max,N_alpha-1);
filepath    = "stock_scraper/combined_history.csv";
days        = 21;

% Options
prnt_last_P_group   = 0;

% Create Portfolio
weights     = zeros(N_stocks, N_alpha);
[weights(:,1), mus, sigmas, P0_vec, names]  = create_portfolio(filepath, alphas(1));

mus     = mus';
sigmas = sigmas';
P0_vec  = P0_vec';
names   = names';
for i = 2:length(alphas)
    [weights(:,i), ~, ~] = create_portfolio(filepath, alphas(i));
end
weights(:,end)  = ones(N_stocks, 1)/N_stocks;
aphas   = alphas.^2;
alphas  = [alphas, 2];

% Get average gainz
gainz_grp   = zeros(N_iter, N_alpha);
gainz       = zeros(N_iter, N_stocks);
% [~, gainz_grp_avg, gainz_avg] = Pb6_iteration(weights, mus ,sigmas, P0_vec, N_alpha, N_stocks, days);
% gainz_grp_min   = gainz_grp_avg;
% gainz_grp_max   = gainz_grp_avg;
for i = 1:N_iter
    [P_group, gainz_grp(i,:), gainz(i,:)] = Pb6_iteration(weights, mus ,sigmas, P0_vec, N_alpha, N_stocks, days);
    % gainz_grp_avg   = gainz_grp_avg + gainz_group;
    % gainz_avg       = gainz_avg + gainz;
end

% Actual average time
gainz_grp_avg   = mean(gainz_grp, 1);
gainz_grp_std   = std(gainz_grp, 1);
gainz_avg       = mean(gainz, 1);
gainz_std       = std(gainz, 1);

% Plotting
if (prnt_last_P_group)
    days_vec    = (1:days+1)';
    figure(1);
    clf;
    plot(days_vec, P_group);
    xlabel("time (days)");
    ylabel("Portfolio Value ($)");
end

% Revenues
figure(2);
clf;
plot(alphas, gainz_grp_avg);
xlabel("Alpha value");
ylabel("Overall Gain ($)");

% Gain Standard Deviations
figure(3);
clf;
plot(alphas, gainz_grp_std);
xlabel("Alpha value");
ylabel("Standard Deviation");

% Analyzing
alpha_low   = alphas(1);
alpha_high  = alphas(end);
[~, ind_max]    = max(gainz_grp_avg);
[~, ind_min]    = min(gainz_grp_avg);
alpha_max   = alphas(ind_max);
alpha_min   = alphas(ind_min);

% Printing
disp("Low Alpha Investment of alpha = " + alpha_low + ": gainz = " + gainz_grp_avg(1));
print_invest_info(names, weights(:,1), gainz_avg, gainz_std, mus, sigmas);
disp("High Risk Investment of alpha = " + alpha_high + ": gainz = " + gainz_grp_avg(end-1));
print_invest_info(names, weights(:,end-1), gainz_avg, gainz_std, mus, sigmas);
disp("Safe Investment of even distribution: gainz = " + gainz_grp_avg(end));
print_invest_info(names, weights(:,end), gainz_avg, gainz_std, mus, sigmas);
disp("Highest gain investment of alpha = " + alpha_max + ": gainz = " + gainz_grp_avg(ind_max));
print_invest_info(names, weights(:,ind_max), gainz_avg, gainz_std, mus, sigmas);
disp("Highest loss investment of alpha = " + alpha_min + ": gainz = " + gainz_grp_avg(ind_min));
print_invest_info(names, weights(:,ind_min), gainz_avg, gainz_std, mus, sigmas);

%% Helpful Functions
function [] = print_invest_info(names, weights, gainz, stds, mus, sigmas)
    % Fixing
    gainz   = gainz';
    stds    = stds';

    % Rounding
    mus     = round(mus,4);
    sigmas  = round(sigmas,4);
    gainz   = round(gainz, 4);
    stds    = round(stds, 4);
    weights     = round(weights, 4);

    % Sorting
    printTable  = table(names, weights, gainz, stds, mus, sigmas);
    printTable  = sortrows(printTable, 3, "descend");

    % Printing
    disp(printTable);
end