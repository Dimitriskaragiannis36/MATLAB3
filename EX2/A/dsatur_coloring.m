function colors = dsatur_coloring(A)
    % Αλγόριθμος DSATUR (Διορθωμένος)
    n = size(A, 1);
    colors = zeros(1, n); 
    degrees = sum(A, 2)'; 
    uncolored = 1:n; 

    while ~isempty(uncolored)
        % 1. Υπολογισμός βαθμού κορεσμού (Πλήθος ΔΙΑΦΟΡΕΤΙΚΩΝ χρωμάτων γειτόνων)
        sat = zeros(1, length(uncolored));
        for i = 1:length(uncolored)
            u = uncolored(i);
            neighbors = find(A(u, :) == 1);
            colored_neighbors = neighbors(colors(neighbors) > 0);
            sat(i) = length(unique(colors(colored_neighbors))); % Μόνο τα unique
        end

        % 2. Εύρεση κόμβων με τον μέγιστο κορεσμό
        max_sat = max(sat);
        candidates = find(sat == max_sat);

        % 3. Tie-breaker: Αν υπάρχουν ισοπαλίες, διάλεξε αυτόν με τον μεγαλύτερο βαθμό
        if length(candidates) > 1
            [~, best_idx] = max(degrees(uncolored(candidates)));
            idx = candidates(best_idx);
        else
            idx = candidates(1);
        end

        node = uncolored(idx);

        % 4. Βρες τα χρώματα που ήδη χρησιμοποιούνται από τους γείτονες
        neighbors = find(A(node, :) == 1);
        used_colors = unique(colors(neighbors(colors(neighbors) > 0)));
        
        % 5. Διάλεξε το μικρότερο διαθέσιμο χρώμα
        color = 1;
        while ismember(color, used_colors)
            color = color + 1;
        end
        
        % Ανάθεση χρώματος
        colors(node) = color;
        
        % Αφαίρεση του κόμβου από τη λίστα των uncolored
        uncolored(idx) = [];
    end
end