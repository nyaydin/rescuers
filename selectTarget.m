function target = selectTarget( nodes, excluded_nodes )

    % TODO: select based on strategy of player.
    
    available_nodes = setdiff(nodes, excluded_nodes);
    target = available_nodes(randi(length(available_nodes)));
end

