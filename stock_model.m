% Models a stock for you
% P0 -- Initial price of the stock
% mu -- Annual Return
% sigma -- Annual Volatility
% delT -- Timestep in years, 1/260 for 1 work day
% days -- How many steps the simulation should take
function [days_vec, P_vec] = stock_model(P0, mu, sigma, delT, days)

    % Attempts to recreate figures such as Fig 12.1 and Fig 12.2 in the book
    % days    = length(days_vec);
    
    % Starting vectors
    P_vec   = ones(1,days+1)*P0;
    days_vec    = 0:days;
    
    % Finding P_vec
    for i = 2:length(P_vec)
        phi     = randn;
        P_t     = P_vec(i-1);
        P_vec(i)    = P_t + mu*P_t*delT + sigma*P_t*sqrt(delT)*phi;
    end

end