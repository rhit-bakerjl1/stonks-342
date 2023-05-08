%% Test Stock Model
close;
clear;
clc;

% Constants
P0  = 100;
sigma   = 0.15;
% sigma = 0.1;
% mu  = 0.03;
mu  = sigma^2/2;
% mu  = 0;
delT    = 1/260;
% days    = 260*1;
days    = 260;

% Options
sim_20_yr   = 0;
plt_example_stocks  = 0;
sim_rand_port = 1;

% Finding average return
if (sim_20_yr) 
    num_sim     = 100000;
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
    N           = 1600;
    sigmaMean   = 0.15;
    % muMean      = 0.03;
    muMean      = 0.03;
    PMean       = 4;
    PMin        = 4;
    xNum        = 15*N;
    [P, x1]     = func_invest_random(N, PMean, PMin, xNum);
    [P3, x3]    = func_invest_random(N, PMean, PMin, xNum);
    [P4, x4]    = func_invest_random(N, PMean, PMin, xNum);
    [P5, x5]    = func_invest_random(N, PMean, PMin, xNum);
    [P6, x6]    = func_invest_random(N, PMean, PMin, xNum);
    [P7, x7]    = func_invest_random(N, PMean, PMin, xNum);
    [sigmas, mus, weights, V]   = func_random_portfol(N, sigmaMean, muMean, P, x1);
    % [sigmas, mus, weights1, V]   = func_struct_portfol(N, sigmaMean, muMean, P, x1);

    % New Investment strategies
    weights2    = ones(size(x1));
    weights2    = weights2/sum(weights2);
    x2          = weights2*V./P;
    
    % Timestepping through Portfolio
    R_mat   = zeros(N, days);
    P_mat   = zeros(N, days);
    P_mat(:,1)  = P;
    t   = 0:delT:(days)*delT;
    for i=1:length(t)-1
        phis    = randn(N,1);
        R_mat(:,i)  = mus*delT + sigmas*sqrt(delT).*phis;
        P_mat(:,i+1)  = P_mat(:,i).*(R_mat(:,i)+1);
    end

    % Overall gains
    P_gain      = P_mat(:,end) - P_mat(:,1);
    P_port1     = x1'*P_mat;
    P_port2     = x2'*P_mat;
    P_port3     = x3'*P_mat;
    P_port4     = x4'*P_mat;
    P_port5     = x5'*P_mat;
    P_port6     = P_mat(1,:)*V/P_mat(1,1);
    P_port7     = P_mat(2,:)*V/P_mat(2,1);

    % Plotting P_mat
    figure(1);
    clf;
    subplot(3,1,1);
    plot(t/delT, P_mat);
    xlabel("Time (days)");
    ylabel("Stock price");
    % legendLabels    = 1:N;
    % legendLabels    = "Stock " + legendLabels;
    % legend(legendLabels);

    % Plotting mus
    subplot(3,1,2);
    plot(mus, P_gain, "o");
    xlabel("Annual Return mu");
    ylabel("Stock Price Gain");

    % Plotting sigmas
    subplot(3,1,3);
    plot(sigmas, P_gain, "o");
    xlabel("Volatility sigma");
    ylabel("Stock Price Gain");

    % Plotting Different Strategies
    figure(2);
    clf;
    plot(t/delT, P_port2);
    hold on;
    plot(t/delT, P_port1);
    plot(t/delT, P_port3);
    plot(t/delT, P_port4);
    plot(t/delT, P_port5);
    plot(t/delT, P_port6);
    plot(t/delT, P_port7);
    xlabel("Time (days)");
    ylabel("Portfolio Value");
    legend("Even Investment", "Random Investment 1", "Random Investment 2", "Random Investment 3", "Random Investment 4", "All in Stock 1", "All in Stock 2");


    % Sample Means, etc. Analyzing overall investment
    P   = P_mat';
    R   = (P(2:end,:)-P(1:end-1,:))./(P(1:end-1,:));
    r   = mean(R);
    C   = cov(R);
    r_vec   = 1/days*sum(R_mat, 2);
    C_mat   = zeros(N, N);
    for i = 1:days
        C_mat   = C_mat + 1/(days-1)*(R_mat(:,i)-mus)*(R_mat(:,i)-mus)';
    end

    % alpha   = (0:0.01:1).^2;
    alpha   = 0.5;
    w   = zeros(N, length(alpha));
    optVal  = zeros(1,length(alpha));
    for i = 1:length(alpha)
        [w(:,i), optVal(:,i)] = quadprog((1-alpha(i))*2*C, -alpha(i)*r, [], [], ...
            ones(1,N), 1, zeros(N,1), ones(N,1));
    end

end

%% Helpful Functions
function [P, x]     = func_invest_random(N, PMean, PMin, xNum)
    % Generate random starting prices
    pSD     = (PMean-PMin)/3;
    P      = pSD*randn(N,1) + PMean;
    % Invest randomly
    xNum    = xNum - N;
    if (xNum < 0) 
        x   = ones(N,1); 
    else
        x   = zeros(N,1);
    end
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
    % sigmas  = linspace(0, sigmaMean*3, sqrt(N))';
    % sigmas  = repmat(sigmas, 1, sqrt(N));
    % sigmas  = reshape(sigmas, N, 1);
    % Generate random mus
    muSD    = muMean/3;
    mus     = muSD*randn(N,1) + muMean;
    % mus     = linspace(0, muMean*3, sqrt(N));
    % mus     = repmat(mus, sqrt(N), 1);
    % mus     = reshape(mus, N, 1);

    % V
    V   = sum(x.*P);
    % Find weights
    weights     = P.*x/V;

end

function [sigmas, mus, weights, V] = func_struct_portfol(N, sigmaMean, muMean, P, x)
    % Generate structured sigmas
    sigmas  = linspace(0, sigmaMean*3, sqrt(N))';
    sigmas  = repmat(sigmas, 1, sqrt(N));
    sigmas  = reshape(sigmas, N, 1);
    % Generate structured mus
    mus     = linspace(0, muMean*3, sqrt(N));
    mus     = repmat(mus, sqrt(N), 1);
    mus     = reshape(mus, N, 1);

    % V
    V   = sum(x.*P);
    % Find weights
    weights     = P.*x/V;
end