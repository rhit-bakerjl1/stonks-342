close;
clear;
clc;

% Constants
P0  = 100;
% sigma   = 0.15;
sigma = 0.1;
% mu  = 0.03;
mu  = 3*sigma^2/2;
% mu  = 0;
delT    = 1/260;
days    = 260*1;

% Options
sim_20_yr   = 1;
plt_example_stocks  = 0;
sim_rand_port = 1;

% Finding average return
if (sim_20_yr) 
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
end

% Plots
if (plt_example_stocks) 
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
end

if (sim_rand_port)
    % Generate random portfolio
    N           = 10;
    sigmaMean   = 0.15;
    muMean      = 0.03;
    PMean       = 4;
    PMin        = 2;
    xNum        = 100;
    [P, x]  = func_invest_random(N, PMean, PMin, xNum);
    [sigmas, mus, weights, V] = func_random_portfol(N, sigmaMean, muMean, P, x);
    
    % Timestepping through Portfolio
    R_mat   = zeros(N, days);
    t   = 0:delT:days/delT;
    for i=1:length(t)
        % R_mat(:,t) = ;
    end

end
    


function [P, x]     = func_invest_random(N, PMean, PMin, xNum)
    % Generate random starting prices
    pSD     = (PMean-PMin)/3;
    P      = pSD*randn(N,1) + PMean;
    % Invest randomly
    x      = ones(N,1);
    xNum    = xNum - N;
    while (xNum > 0)
        addTo   = randi([1,N]);
        xNum = xNum - 1;
        x(addTo)   = x(addTo) + 1;
    end
end

function [sigmas, mus, weights, V] = func_random_portfol(N, sigmaMean, muMean, P, x)
    % Generate random sigmas
    sigmaSD     = sigmaMean/3;
    sigmas  = sigmaSD*randn(N,1) + sigmaMean;
    % Generate random mus
    muSD    = muMean/3;
    mus     = muSD*randn(N,1) + muMean;

    % V
    V   = sum(x.*P);
    % Find weights
    weights     = P.*x/V;

end
