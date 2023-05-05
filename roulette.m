clear;
close;
clc;

n = 1000000;
start_money = 2000;
money = zeros(n, 1);

% For each simulation
for sim = 1 : n
    % For each year
    m = start_money;
    for month = 1 : 6
        % Spin the wheel
        ball = randi([1 38]);
        if ball <= 18
            m = m * 2;
        else
            m = 0;
        end
    end
    money(sim) = m;
end

mean(money)
max(money)
std(money)
groupcounts(money)