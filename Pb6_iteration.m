function [P_group, P_mat, gainz_group, gainz] = Pb6_iteration(weights, mus, sigmas, alphas, P0_vec, N_alpha, N_stocks, days, daisy)
    % Constants
    delT        = 1/260;

    % Process my weights
    V       = 20000;     % Dollars we invest
    P0_mat  = repmat(P0_vec, 1, N_alpha);
    
    % Simulate Stocks with daisychaining
    for iter = 1:4
        if (iter ~= 1)
            [weights(:,1), mus, sigmas, P0_vec]  = create_portfolio_p(P_mat, alphas(1), days);
            for i = 2:length(alphas)
                [weights(:,i), ~, ~, ~] = create_portfolio_p(P_mat, alphas(i), days);
            end
            weights(:,end)  = ones(N_stocks, 1)/N_stocks;
            mus     = mus';
            sigmas = sigmas';
            P0_vec  = P0_vec';
        end
        phis    = randn(N_stocks, days);
        mus_mat = repmat(mus,1,days);
        sigmas_mat  = repmat(sigmas,1,days);
        R_mat   = mus_mat*delT + sigmas_mat*sqrt(delT).*phis;
        
        % Find how much we got oh boy yeah less go rakin in the dough come on now
        P_group     = zeros(days+1,N_alpha);
        P_mat       = zeros(days+1,N_stocks);
        P_mat(1,:)  = P0_vec;
        P_group(1,:)  = V;
        R_group     = (weights'*R_mat)';
        for day = 1:days
           P_mat(day+1, :)      = P_mat(day, :).*(R_mat(:, day)' + 1);
           P_group(day+1, :)    = P_group(day, :).*(R_group(day, :)+1); 
        end
    
        % Revenues
        gainz   = P_mat(end,:) - P_mat(1,:);
        gainz_group   = P_group(end,:) - P_group(1,:);

        % Updating V
        V   = V + gainz_group;

        if (~daisy)
            break;
        end
    end
end