%% Analyzing Portfolio Choices
close;
clear;
clc;

% Setup for create_portfolio
N_alpha     = 100;
N_stocks    = 44;
% N_iter      = 100000;
N_iter      = 20;
alpha_min   = 0;
alpha_max   = 1;
alphas      = linspace(alpha_min,alpha_max,N_alpha-1);
filepath    = "stock_scraper/combined_history.csv";
days        = 21;

% Options
prnt_last_P_group   = 0;
daisy   = 1;

% Create Portfolio
weights     = zeros(N_stocks, N_alpha);
alphas   = alphas.^2;
alphas  = [alphas, 2];
[weights(:,1), mus, sigmas, P0_vec, names]  = create_portfolio(filepath, alphas(1));

for i = 2:length(alphas)
    [weights(:,i), ~, ~] = create_portfolio(filepath, alphas(i));
end
% New Strategy: invest equally
weights(:,end)  = ones(N_stocks, 1)/N_stocks;

% Flip because we said so
mus     = mus';
sigmas = sigmas';
P0_vec  = P0_vec';
names   = names';

% Get average gainz
gainz_grp   = zeros(N_iter, N_alpha);
gainz       = zeros(N_iter, N_stocks);
% [~, gainz_grp_avg, gainz_avg] = Pb6_iteration(weights, mus ,sigmas, P0_vec, N_alpha, N_stocks, days);
% gainz_grp_min   = gainz_grp_avg;
% gainz_grp_max   = gainz_grp_avg;
for i = 1:N_iter
    [P_group, P_mat, gainz_grp(i,:), gainz(i,:)] = Pb6_iteration(weights, mus ,sigmas, alphas, P0_vec, N_alpha, N_stocks, days, daisy);
    % gainz_grp_avg   = gainz_grp_avg + gainz_group;
    % gainz_avg       = gainz_avg + gainz;
    disp("i = " + i);
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
plot(alphas(1:end-1), gainz_grp_avg(1:end-1), "LineWidth", 2.5);
xlabel("Risk Level (alpha)");
ylabel("Average Profit after One Month ($)");

% Gain Standard Deviations
figure(3);
clf;
errorbar(alphas(1:end-1), gainz_grp_avg(1:end-1) + 2000, gainz_grp_std(1:end-1));
ylim([0 3000]);
xlabel("Risk Level (alpha)");
ylabel("Portfolio Value after One Month ($)");

% Finding percentage chance of going positive
zero_std    = -gainz_grp_avg./gainz_grp_std;
pos_percents    = zeros(size(zero_std));
for i = 1:N_alpha
    pos_percents(i) = diff(normcdf([zero_std(i), 50]));
end

% Plotting percentage chance of going positive
figure(4);
clf;
plot(alphas, pos_percents*100, "LineWidth", 2.5);
xlabel("Alpha value");
ylabel("Chance of going positive (%)");
xlim([0,1]);

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