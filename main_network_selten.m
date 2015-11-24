% There is a normal network, which gets perturbed by a calamity.
% Rescuers explores the new perturbed network trying to rescue people.
% Rescuers start from hospital (node 1) and must find and bring back 
% injured people before it is too late.
% Different strategies of exploration and rescuing are tested.
% Goal: find out if there is a superior strategy or a mix thereof.


% TODO: Example of strategies for rescuers: (no need to implement all)

% Strategy 0: use info on perturbed network. Baseline.
% Strategy 1: use information only from known network (not the perturbed)
% Strategy 2: update known information of network during exploration
% Strategy 3: exchange information with other players.
% Strategy 4: explore all network first then start rescuing people
% How is the payoff of a rescuer computed?
% Get points if an injured is found.
% How are the propensity adjusted?

clc
close all
clear

% Open debugger if there is an error.
% dbstop if error


%% Create Networks

from_nodes = [6 1 2 2 3 4 4 5 5 6 1];
to_nodes = [2 6 3 5 4 1 6 3 4 3 5];

% Use discrete weights, so we can consider them a single simulation steps.
Wt = [4 10 1 12 14 8 8 9 7 7 9]; %weights of time
Nt = sparse(from_nodes, to_nodes, Wt); %sparse matrix that is the network

% view(biograph(Nt,[],'ShowWeights','on')) %display the network with weights

% Perturn the network.
Ws = Wt * (1+rand()); 
Ns = sparse(from_nodes, to_nodes, Ws); %sparse matrix that is the network

% view(biograph(Ns,[],'ShowWeights','on')) %display the network with weights

% graphallshortestpaths(Ns) %evaluate total time spent from every node of the network to each other node

% Create biograph object (used to compute shortest paths).
Obj_s = biograph(Ns); 
Obj_t = biograph(Nt); 

% Nr. Rescuers
N = 10;

% Nr. Rounds = the time to live (ttl). Injured people need to be found,
% and brought to hospital in less than this time, otherwise they die.
% This is a simplication. All injured people have the same time to live.
T = 200;

% Nr. of rescued people in the simulation.
RESCUED = 0;
TOTAL_INJURED = 0; % will be init later.

% Multiplier for how many injured per node maximum.
INJ_MULTIPLIER = 20;

% Under normal conditions, the time to reach a node from another one.
normal_sp = graphallshortestpaths(Nt);

% Under a natural disaster the network was perturbed, and rescuers do
% not know for sure how much it takes to move between nodes.
% The values of this matrix have to be learnt by the agents (rescuers).
perturbed_sp = graphallshortestpaths(Ns);

% Each rescuer build an own image of the network after perturbation.
% It starts as a copy of the normal network, equal for all rescuers.
% It might get updated during the exploration or not, depending on
% strategy.
networks = repmat(normal_sp, N);

% Number of nodes in the network.
nNodes = length(Nt);

% Node indexes.
nodes = 1:nNodes;

% Place injured people randomly on the network.
injured = (randn(nNodes,1) < 0.5) .* randi(INJ_MULTIPLIER, nNodes, 1);

% Place Hospital in the network (node 1);
injured(1) = 0;

TOTAL_INJURED = sum(injured);

% Current position in the network for the rescuers.
% (At the beginning all in position 1 - hospital).
positions = ones(N,1);

% Reformat network data for use in the main simulation loop.
% If a rescuer is in node X and wants to reach node Y, which node should
% he / she go to next, accordint to the shortest path algorithm?
next_move = zeros(nNodes);
for x = 1 : nNodes
    for y = 1 : nNodes
        if (x ~= y) 
            [dist_t, path_t, pred_t] = shortestpath(Obj_t, x,y);
            next_move(x,y) = path_t(2);
        end
    end
end

% Where a rescuer is heading to. 0 means no target.
target_nodes = zeros(N, 1);

% Whether a rescuer is carrying an injured.
carrying_injured = zeros(N, 1);

% Current waiting time in node.
% To move from one node to another, rescuers need to wait the value
% of the wait in the connection between two nodes.
waiting_in_node = zeros(N,1);

% Available strategies.
avail_strategies = [ 1 2 3 4 ];

% Nr. available strategies each player each round.
nr_strategies = length(avail_strategies);

% Currently chosen strategy by each player.
strategies = ones(N,1);

% Initial propensities for each strategy of each player.
propensities = ones(N, nr_strategies);

% Probabilities to choose a strategy. (all equally probable at t=0).
probabilities = ones(N, nr_strategies)*0.25;

% Records payoffs and routes for each player at each iteration.
payoff_matrix = zeros(N,T);
route_choice = zeros(N,T);

%% Model the strategy choices of each player 
% For each round, each player decides which strategy to choose
% depending on the previous payoff and propensity. 
% Probability is given by x(i,j)/sum(x(i,j))

    

%% Routh choice and interaction
for t = 1 : T
        
        for player = 1 : N
            
            % Reset variables for each player.
            injured_found = 0;
            injured_rescued = 0;
            
            % Randomly draw a strategy according to the prob distribution.
            curStrategy = randsample(avail_strategies, ...
                1, true, probabilities(player,:));
            
            strategies(player) = curStrategy;
            
            % Find a target node (if none was assigned before).
            % Each rescuer randomly decide to go to a node and check
            % if there are people to rescue there.
            % TODO: this could be decided according to the strategy
            % currently played by player.
            if (target_nodes(player) == 0)
                % Notice: +1 is to exclude the hospital as a target.
                target_nodes(player) = selectTarget(nodes, ...
                                                    [positions(player)]);
                                                
                [next_node, time_to_node] = enterNewNode(player, ...
                                                         positions, ...
                                                         next_move, ...
                                                         target_nodes, ...             
                                                         perturbed_sp);
                                                     
                                                     
                % Set waiting time in current node.
                waiting_in_node(player) = time_to_node;
                                                     
            end
            
            % Entering a new node.
            if (waiting_in_node(player) <= 0)
                
                [next_node, time_to_node] = enterNewNode(player, ...
                                                         positions, ...
                                                         next_move, ...
                                                         target_nodes, ...             
                                                         perturbed_sp);

                % Move player in new position.
                positions(player) = next_node;
                
                % Player is in hospital
                if (next_node == 1)
                    
                    % Transit or rescued successful?
                    if (carrying_injured(player) == 1)
                        injured_rescued = 1;
                        carrying_injured(player) = 0;
                        RESCUED = RESCUED + 1;
                        % TODO: you might want to save the time as well.
                        % for more statistics.
                        % TODO: you might want to save the strategy type
                        % who saved the injured.
                        
                        % Reset target node and time_to_node.
                        % They will be initialized at next iteration.
                        target_nodes(player) = 0;
                        time_to_node = 0;
                    end
                end
                
                % Evaluate node.
                
                % TODO: some strategies might not look for injured
                % people in the middle of their travel, i.e. this check
                % is performed only if next_node = target_node.
                % Limits to 1 person rescued at the same time.
                if (injured(next_node) > 0 && ...
                        carrying_injured(player) == 0)
                    
                    % If an injured is found, remove it,
                    % change target_node to hospital,
                    % and update other variables.
                    injured(next_node) = injured(next_node) - 1;
                    carrying_injured(player) = 1;
                    target_nodes(player) = 1;
                    injured_found = true;
                    [next_node, time_to_node] = enterNewNode(player, ...
                                                         positions, ...
                                                         next_move, ...
                                                         target_nodes, ...             
                                                         perturbed_sp);
                else
                    if (next_node == target_nodes(player))
                        % We reached target node, but did not find an
                        % injured person. Reset!
                        target_nodes(player) = 0;
                        time_to_node = 0;
                    end
                end
                
                % Set waiting time in current node.
                waiting_in_node(player) = time_to_node;
            else
                % Reduce waiting time.
                waiting_in_node(player) = waiting_in_node(player) - 1;    
            end
            
          
            % Compute payoff each player.
            % TODO: Need to decide how to do it exactly. 
            % Might be strategy dependent or not.
            payoff = 0;
            if (injured_found)
                payoff = payoff + 50;
            elseif (injured_rescued)
                payoff = payoff + 100;
            else
                payoff = payoff - 2;
            end
            
            % Updates propensities.
            if payoff >= 0
                propensities(curStrategy) = ...
                    propensities(curStrategy) + payoff;
            else
                not_choosen_strategies = ...
                    setdiff(avail_strategies, curStrategy);
                propensities(not_choosen_strategies) = ...
                    propensities(not_choosen_strategies) - payoff;
            end
            
        end
             
end


RESCUED
DEAD = TOTAL_INJURED - RESCUED


