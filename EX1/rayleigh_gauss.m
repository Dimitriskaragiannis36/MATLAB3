function [lambda, znew, iter] = rayleigh_gauss(A, y, q, tol, maxiter)
    % Μέθοδος πηλίκων Rayleigh με αρχική μετατόπιση q
    n = length(A);
    I = eye(n);
    lambda0 = q; 
    
    % Κανονικοποίηση αρχικού διανύσματος με 2-νόρμα (απαραίτητο για Rayleigh)
    z = y / norm(y); 
    
    for iter = 1:maxiter
        % Επίλυση του (A - qI)w = z
        w = gauss_elimination(A - q * I, z);
        
        % Κανονικοποίηση νέου διανύσματος
        znew = w / norm(w);
        
        % Υπολογισμός ιδιοτιμής μέσω του πηλίκου Rayleigh
        lambda = znew' * A * znew;
        
        % Έλεγχος σύγκλισης
        if abs(lambda - lambda0) < tol
            return;
        end
        
        % Ενημέρωση του q (Αυτή είναι η ουσία της μεθόδου Rayleigh)
        q = lambda;
        lambda0 = lambda;
        z = znew;
    end
    warning('Όχι σύγκλιση μετά από %d επαναλήψεις', maxiter);
end