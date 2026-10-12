function x = gauss_elimination(A, b)
    n = length(b);
    a = [A, b]; 
    
    % Απαλοιφή με μερική οδήγηση
    for k = 1:n-1
        [~, max_idx] = max(abs(a(k:n, k)));
        max_idx = max_idx + k - 1;
        
        if max_idx ~= k
            a([k max_idx], :) = a([max_idx k], :);
        end
        
        for i = k+1:n
            m = a(i, k) / a(k, k);
            a(i, k:end) = a(i, k:end) - m * a(k, k:end);
        end
    end
    
    % Οπισθοδιαγραφή
    x = zeros(n, 1);
    x(n) = a(n, end) / a(n, n);
    for i = n-1:-1:1
        x(i) = (a(i, end) - a(i, i+1:n) * x(i+1:n)) / a(i, i);
    end
end