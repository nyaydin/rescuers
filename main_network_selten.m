% Main_network_selten: 
% a network is created with two different sets of weights (the time weights
% and the distance weights)--> the "time network" and the "distance network". The shortest path is evaluated for the time
% network and the distance network. Then the Selten agent based model is
% run, the stategies are: 1) choose the shortest time path to reach a node, 2)
% choose the shortest distance path to reach a node, 3) direct strategy,
% 4)contrarian strategy.
clc
close all
clear


%% Create a network of time
Wt = [.41 .99 .51 .32 .15 .88 .38 .32 .36 .29 .21]; %weights of time
Nt= sparse([6 1 2 2 3 4 4 5 5 6 1],[2 6 3 5 4 1 6 3 4 3 5],Wt) %sparse matrix that is the network

view(biograph(Nt,[],'ShowWeights','on')) %display the network with weights

graphallshortestpaths(Nt) %evaluate total time spent from every node of the network to each other node

Obj_t=biograph(Nt); %in order to evaluate the shortest path
[dist_t, path_t, pred_t] = shortestpath(Obj_t, 1,2); %in path the node of the shortest path are displayed


%Create a network of distance
Ws = [.16 .18 .12 .20 .15 .11 .12 .13 .14 .15 .16]; %weights of distance
Ns = sparse([6 1 2 2 3 4 4 5 5 6 1],[2 6 3 5 4 1 6 3 4 3 5],Ws) %sparse matrix that is the network

view(biograph(Ns,[],'ShowWeights','on')) %display the network with weights

graphallshortestpaths(Ns) %evaluate total time spent from every node of the network to each other node

Obj_s=biograph(Ns); %in order to evaluate the shortest path
[dist_s, path_s, pred_s] = shortestpath(Obj_s, 1,2); %in path the node of the shortest path are displayed

%% Evaluate with Selten model 

% Main road: shortest path of the time network
% Secondary road: shortest path of the distance network

%% Define main variables
% Nr. Players --> Rescuers
N = 18;
% Nr. Rounds.
T = 200;
% Initial endowment of each player.
A0 = 200;
% Initial endowment vector for all players.
A0_vec = A0*ones(N,1);
% Initial propensities for each strategy of each player.
propens_vec_0 = [1 2];
% Nr. available strategies each player each round.
nr_strategies = 2;
% Vector position strategies.
position_strategy_vec = 1 : nr_strategies;
% Max payoff each player
max_payoff = 40;
% Penalty for time spent on main road (M).  
Pen_M = dist_t;
% Penalty for time spent on secondary road (S).
Pen_S = dist_s;
% Multiplier penalty for time spent on main road.
Multi_M = 2;
% Multiplier penalty for time spent on secondary road.
Multi_S = 3;


%% Model the strategy choices of each player 
% For each round, each player decides which strategy to choose
    % depending on the previous payoff and propensity. 
    % Probability is given by x(i,j)/sum(x(i,j))

% Variables at t1
strategy_t1 = zeros(N,1);
payoff_vector_t1 = zeros(N,1);
propensity_vector_before_draw_t1 = zeros(N,1);

% At initial time only strategies 1 (Main) and 2 (Secondary) are available
strategies_vec_t1 = propens_vec_0(1:2);
probability_strategies_vec_t1 = strategies_vec_t1/sum(strategies_vec_t1);
position_strategy_vec_t1 = position_strategy_vec(1:2);

% Variables for t>1
propensity_array = zeros(T,nr_strategies,N);
propensity_array_0 = repmat(propens_vec_0,1,1,N); 
propensity_array(1,:,:) = propensity_array_0;
propensity_cell = cell(1,N);
payoff_matrix = zeros(N,T);
route_choice = zeros(N,T); % CHANGE THIS LATER, pay attention on this

%% Routh choice and interaction
for t = 1 : T 
    if t == 1
        
        for player = 1 : N          
            % Randomly draw a strategy (position) according
                %to the probability distribution
            strategy_t1(player) = randsample(position_strategy_vec_t1,1,true,probability_strategies_vec_t1);
        end

        %%%%%%%% Payoff of each player depends on  
        % how many players have chosen each strategy
        % Nr. players choosing strategy 1 (Main)
        n_M_t1 = sum(strategy_t1 == 1);
        % Nr. players choosing strategy 2 (Secondary)
        n_S_t1 = sum(strategy_t1 == 2);
        % Penalty for players choosing strategy 1 (Main)
        tM_t1 = Pen_M + Multi_M*n_M_t1;
        % Penalty for players choosing strategy 2 (Secondary)
        tS_t1 = Pen_S + Multi_S*n_S_t1;
        % Payoff for players choosing strategy 1 (Main)
        a_M_t1 = max_payoff - tM_t1;
        % Payoff for players choosing strategy 2 (Secondary)
        a_S_t1 = max_payoff - tS_t1;
        
        % Payoff vector of the 18 players
        payoff_vector_t1(strategy_t1 == 1) = a_M_t1;
        payoff_vector_t1(strategy_t1 == 2) = a_S_t1;
        
        % Initial propensity of each player according to draw at t1
        propensity_vector_before_draw_t1(strategy_t1 == 1) = strategies_vec_t1(1);
        propensity_vector_before_draw_t1(strategy_t1 == 2) = strategies_vec_t1(2);
        % Compute the new propensities of each player
        for player = 1 : N
            payoff = payoff_vector_t1(player);
            strategy = strategy_t1(player);
            propensity_mat = propensity_array(:,:,player);
            not_choosen_strategies = setdiff(position_strategy_vec,strategy);
            if payoff >= 0
                propensity_mat(1,strategy) = propensity_mat(1,strategy) + payoff;
            else
                propensity_mat(1,not_choosen_strategies) = ...
                    propensity_mat(1,not_choosen_strategies) - payoff;
            end
            % Replace with the new propensities
            propensity_array(1,:,player) = propensity_mat(1,:);       
            % Store route choice of each player for further rounds
            route_choice(:,1) = strategy_t1;
            
            % Save payoffs each player
            payoff_matrix(:,1) = payoff_vector_t1;
            
        end      
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For the remaining time periods          
    else
        for player = 1 : N
            % Recompute the probabilities to choose each strategy 
            % based on previous step
            all_previous_propensities = propensity_array(:,:,player);
            propensity_t = all_previous_propensities(t-1,:);
            route_choice_t = route_choice(player,t-1);           
            probability_strategies_vec = propensity_t/sum(propensity_t);
            payoff_t = payoff_matrix(player,t-1);
            all_prev_payoffs = payoff_matrix(player,1:t-1);
            median_all_prev_payoffs = median(all_prev_payoffs);
            count_smaller_median = sum(median_all_prev_payoffs > all_prev_payoffs);
            count_bigger_median = sum(median_all_prev_payoffs < all_prev_payoffs);
            
            % Choose new strategy according to new propensity
            strategy_tp1 = randsample(position_strategy_vec,1,...
                                true,probability_strategies_vec);
            strategy(player) = strategy_tp1;
            % If new strategy is 1 (Main)
            if strategy_tp1 == 1
                route_choice(player,t) = strategy_tp1;
            % If the new strategy is 2 (Secondary)
            elseif strategy_tp1 == 2
                route_choice(player,t) = strategy_tp1;
            elseif strategy_tp1 == 3
                % If the payoff is lower than before, change route
                if payoff_t < median_all_prev_payoffs
                    if route_choice_t == 1
                        route_choice(player,t) = 2;
                    else
                        route_choice(player,t) = 1;
                    end
                % If the payoff is higher than before, keep route
                elseif payoff_t > median_all_prev_payoffs
                    route_choice(player,t) = route_choice_t;
                % If payoff is the same as before, decision rule depends on
                    % nr. of occurances
                else
                    if count_smaller_median < count_bigger_median
                        route_choice(player,t) = route_choice_t;
                    elseif count_smaller_median > count_bigger_median
                        if route_choice_t == 1
                            route_choice(player,t) = 2;
                        else
                            route_choice(player,t) = 1;
                        end
                    else
                        route_choice(player,t) = randsample([1 2],1,true,[0.5 0.5]);
                    end
                end
            % If selected strategy is nr. 4 (Contrarian)
            else
                if payoff_t > median_all_prev_payoffs
                    if route_choice_t == 1
                        route_choice(player,t) = 2;
                    else
                        route_choice(player,t) = 1;
                    end
                elseif payoff_t < median_all_prev_payoffs
                    route_choice(player,t) = route_choice_t;
                else
                    if count_smaller_median > count_bigger_median
                        route_choice(player,t) = route_choice_t;
                    elseif count_smaller_median < count_bigger_median
                        if route_choice_t == 1
                            route_choice(player,t) = 2;
                        else
                            route_choice(player,t) = 1;
                        end
                    else
                        route_choice(player,t) = randsample([1 2],1,true,[0.5 0.5]);
                    end
                end  
            end
        end
        
        % Compute payoff each player
        % Nr. players choosing strategy 1 (Main)
        n_M = sum(route_choice(:,t) == 1);
        % Nr. players choosing strategy 2 (Secondary)
        n_S = sum(route_choice(:,t) == 2);
        % Penalty for players choosing strategy 1 (Main)
        tM = Pen_M + Multi_M*n_M;
        % Penalty for players choosing strategy 2 (Secondary)
        tS = Pen_S + Multi_S*n_S;
        % Payoff for players choosing strategy 1 (Main)
        a_M = max_payoff - tM;
        % Payoff for players choosing strategy 2 (Secondary)
        a_S = max_payoff - tS;
        
                
        % Payoff vector of the 18 players
        payoff_matrix(route_choice(:,t) == 1,t) = a_M;
        payoff_matrix(route_choice(:,t) == 2,t) = a_S;
        
        % Update propensity
        for player = 1 : N
            payoff = payoff_matrix(player,t);
            strategy_now = strategy(player);
            propensity_mat = propensity_array(:,:,player);
            propensity_mat(t,:) = propensity_mat(t-1,:);
            not_choosen_strategies = setdiff(position_strategy_vec,strategy_now);
            if payoff >= 0
                propensity_mat(t,strategy_now) = propensity_mat(t-1,strategy_now) + payoff;
            else
                propensity_mat(t,not_choosen_strategies) = ...
                    propensity_mat(t-1,not_choosen_strategies) - payoff;
            end
            % Replace with the new propensities
            propensity_array(t,:,player) = propensity_mat(t,:);
        end
    end
            
             
end


%% Plots

% 1. Number of players choosing strategy 1 or 2 each round
nr_players_main = sum(route_choice == 1);
max_nr_players_main = max(nr_players_main);
vec_max_nr_players_main = max_nr_players_main*ones(1,T);
min_nr_players_main = min(nr_players_main);
vec_min_nr_players_main = min_nr_players_main*ones(1,T);
ave_nr_players_main = mean(nr_players_main);
vec_ave_nr_players_main = ave_nr_players_main*ones(1,T);

nr_players_secondary = sum(route_choice == 2);
max_nr_players_secondary = max(nr_players_secondary);
vec_max_nr_players_secondary = max_nr_players_secondary*ones(1,T);
min_nr_players_secondary = min(nr_players_secondary);
vec_min_nr_players_secondary = min_nr_players_secondary*ones(1,T);
ave_nr_players_secondary = mean(nr_players_secondary);
vec_ave_nr_players_secondary = ave_nr_players_secondary*ones(1,T);

time_vec = 1 : T;

ylim_plot = [-2 20];
xlim_plot = [1 T];

figure
subplot(2,1,1)
plot(time_vec,nr_players_main)
hold on
plot(time_vec,[vec_max_nr_players_main;...
               vec_min_nr_players_main;...
               vec_ave_nr_players_main])
ylim(ylim_plot)
xlim(xlim_plot)
title('Number of players choosing the main road')
legend('Nr. players','Max nr. players','Min nr. players','Average nr. players')
xlabel('Round')

subplot(2,1,2)
plot(time_vec,nr_players_secondary)
title('Number of players choosing the secondary road')
hold on
plot(time_vec,[vec_max_nr_players_secondary;...
               vec_min_nr_players_secondary;...
               vec_ave_nr_players_secondary])
ylim(ylim_plot)
xlim(xlim_plot)
xlabel('Round')
legend('Nr. players','Max nr. players','Min nr. players','Average nr. players')


%% Scatter plot nr. road changes vs. cumulative payoff
cum_payoff = cumsum(payoff_matrix,2);
sum_cum_payoff = cum_payoff(:,end);

route_changes = zeros(N,T-1);
for t = 2 : T
    route_choice_t = route_choice(:,t);
    route_choice_tm1 = route_choice(:,t-1);
    route_change_indie_mat = route_choice_t ~= route_choice_tm1;
    route_changes(:,t-1) = route_change_indie_mat;
end

route_changes_total = sum(route_changes,2);

figure
% Plot points.
scatter(route_changes_total,sum_cum_payoff)
xlabel('Nr. route changes each player')
ylabel('Cumulative payoffs each player')
title('Link between route changes and cumulative payoff')
hold on;
% Plot regression line.
my_poly=polyfit(route_changes_total,sum_cum_payoff,1);
X2 = 1:max(route_changes_total+10); % X data range
Y2 = polyval(my_poly,X2);
plot(X2,Y2);









