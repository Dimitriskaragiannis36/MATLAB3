% Ορισμός των πινάκων A1 και A2 (Διορθώθηκε η συμμετρία του A1)
A1 = [0 1 1 1 1 0 0 0; 
      1 0 1 0 0 0 1 0;
      1 1 0 1 0 1 0 0;
      1 0 1 0 1 1 0 0;
      1 0 0 1 0 0 0 0;
      0 0 1 1 0 0 1 1;
      0 1 0 0 0 1 0 1; % <-- Διόρθωση συμμετρίας
      0 0 0 0 0 1 1 0];

A2 = [0 1 1 0 0 0 0 0 0; 
      1 0 1 0 0 0 0 1 0;
      1 1 0 0 1 0 0 0 1;
      0 0 0 0 1 1 1 0 0;
      0 0 1 1 0 0 1 0 0;
      0 0 0 1 0 0 1 0 0;
      0 0 0 1 1 1 0 0 1;
      0 1 0 0 0 0 0 0 1;
      0 0 1 0 0 0 0 1 0];

matrices = {A1, A2};
matrix_names = {'A1', 'A2'};

% Παράμετροι
tol = 1/2 * 10^-6;
maxiter = 100;

for m = 1:2
    A = matrices{m};
    n = size(A, 1);
    y = ones(n, 1);
    
    disp(' ');
    disp(['========================================================']);
    disp(['                 ΑΝΑΛΥΣΗ ΓΡΑΦΗΜΑΤΟΣ ', matrix_names{m}]);
    disp(['========================================================']);

    % Εύρεση ιδιοτιμών και σωστή ταξινόμηση κατά ΑΠΟΛΥΤΗ ΤΙΜΗ
    D = eig(A);
    [~, idx] = sort(abs(D), 'descend');
    sorted_eig = D(idx);
    
    l1 = sorted_eig(1); % Μεγαλύτερη κατά μέτρο
    l2 = sorted_eig(2); % 2η μεγαλύτερη κατά μέτρο
    ln = sorted_eig(end); % Μικρότερη κατά μέτρο

    fprintf('Μεγαλύτερη ιδιοτιμή (λ1): %.4f\n', l1);
    fprintf('Δεύτερη μεγαλύτερη (λ2): %.4f\n', l2);
    fprintf('Μικρότερη ιδιοτιμή (λn): %.4f\n\n', ln);

    % ---- ΥΠΟΛΟΓΙΣΜΟΣ λ1 ----
    if m == 1
        q_values1 = 2.5:0.1:3.7;
    else
        q_values1 = 1.8:0.1:3.0;
    end
    
    disp('--- Υπολογισμός λ1 (Αντίστροφη Μέθοδος Δυνάμεων) ---');
    fprintf('%-10s | %-12s | %-10s | %-20s\n', 'Μετατόπιση', 'Επαναλήψεις', 'Ιδιοτιμή', 'Ιδιοδιάνυσμα');
    t_power_l1 = tic;
    for q = q_values1
        [lambda_approx, z, itcount] = eig_power_gauss(A, y, q, tol, maxiter);
        z_str = sprintf('[%.2f, %.2f, %.2f, %.2f...]', z(1), z(2), z(3), z(4));
        fprintf('%-10.2f | %-12d | %-10.4f | %s\n', q, itcount, lambda_approx, z_str);
    end
    t_power_l1 = toc(t_power_l1);

    disp('--- Υπολογισμός λ1 (Μέθοδος Rayleigh) ---');
    fprintf('%-10s | %-12s | %-10s | %-20s\n', 'Μετατόπιση', 'Επαναλήψεις', 'Ιδιοτιμή', 'Ιδιοδιάνυσμα');
    t_ray_l1 = tic;
    for q = q_values1
        [lambda_approx, z, itcount] = rayleigh_gauss(A, y, q, tol, maxiter);
        z_str = sprintf('[%.2f, %.2f, %.2f, %.2f...]', z(1), z(2), z(3), z(4));
        fprintf('%-10.2f | %-12d | %-10.4f | %s\n', q, itcount, lambda_approx, z_str);
    end
    t_ray_l1 = toc(t_ray_l1);

    % ---- ΥΠΟΛΟΓΙΣΜΟΣ λn ----
    if m == 1
        q_valuesn = -1.1:-0.1:-1.4;
    else
        q_valuesn = -1.3:0.6:-0.8;
    end
    
    disp(' ');
    disp('--- Υπολογισμός λn (Αντίστροφη Μέθοδος Δυνάμεων) ---');
    fprintf('%-10s | %-12s | %-10s | %-20s\n', 'Μετατόπιση', 'Επαναλήψεις', 'Ιδιοτιμή', 'Ιδιοδιάνυσμα');
    t_power_ln = tic;
    for q = q_valuesn
        [lambda_approx, z, itcount] = eig_power_gauss(A, y, q, tol, maxiter);
        z_str = sprintf('[%.2f, %.2f, %.2f, %.2f...]', z(1), z(2), z(3), z(4));
        fprintf('%-10.2f | %-12d | %-10.4f | %s\n', q, itcount, lambda_approx, z_str);
    end
    t_power_ln = toc(t_power_ln);

    disp('--- Υπολογισμός λn (Μέθοδος Rayleigh) ---');
    fprintf('%-10s | %-12s | %-10s | %-20s\n', 'Μετατόπιση', 'Επαναλήψεις', 'Ιδιοτιμή', 'Ιδιοδιάνυσμα');
    t_ray_ln = tic;
    for q = q_valuesn
        [lambda_approx, z, itcount] = rayleigh_gauss(A, y, q, tol, maxiter);
        z_str = sprintf('[%.2f, %.2f, %.2f, %.2f...]', z(1), z(2), z(3), z(4));
        fprintf('%-10.2f | %-12d | %-10.4f | %s\n', q, itcount, lambda_approx, z_str);
    end
    t_ray_ln = toc(t_ray_ln);

    % ---- ΣΥΝΟΨΗ ΧΡΟΝΩΝ (Διορθωμένα Ονόματα Μεταβλητών) ----
    disp(' ');
    fprintf('Χρόνος eig_power (λ1): %.6f s | rayleigh (λ1): %.6f s\n', t_power_l1, t_ray_l1);
    fprintf('Χρόνος eig_power (λn): %.6f s | rayleigh (λn): %.6f s\n', t_power_ln, t_ray_ln);
    disp('--------------------------------------------------------');
end