function [lambda, znew, iter] = eig_power_gauss(A, y, q, tol, maxiter)
    % Τροποποιημένη Αντίστροφη Μέθοδος Δυνάμεων
    n = length(A);
    I = eye(n);
    
    % Αρχική κανονικοποίηση με νόρμα απείρου
    [~, k] = max(abs(y));
    z = y / y(k); 
    
    lambda0 = 0;

    for iter = 1:maxiter
        % Επίλυση του (A - qI)w = z 
        % Χρησιμοποιούμε τη συνάρτηση gauss_elimination που φτιάξατε (ή την backslash \ του MATLAB/Octave)
        w = gauss_elimination(A - q * I, z);
        
        % Βρίσκουμε το στοιχείο του w με τη μέγιστη απόλυτη τιμή
        [~, k] = max(abs(w));
        mu = w(k);
        
        % Υπολογισμός νέου διανύσματος
        znew = w / mu;
        
        % Υπολογισμός προσεγγιστικής ιδιοτιμής
        lambda = q + 1 / mu;
        
        % Έλεγχος σύγκλισης
        if abs(lambda - lambda0) < tol
            return;
        end
        
        z = znew;
        lambda0 = lambda;
    end
    warning('Όχι σύγκλιση μετά από %d επαναλήψεις', maxiter);
end