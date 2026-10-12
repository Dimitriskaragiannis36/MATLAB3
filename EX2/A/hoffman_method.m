% Ορισμός του πίνακα γειτνίασης για το 1ο γράφημα (Διορθώθηκε η συμμετρία)
A1 = [0 1 1 1 1 0 0 0;
      1 0 1 0 0 0 1 0;
      1 1 0 1 0 1 0 0;
      1 0 1 0 1 1 0 0;
      1 0 0 1 0 0 0 0;
      0 0 1 1 0 0 1 1;
      0 1 0 0 0 1 0 1; % <--- Διορθώθηκε το στοιχείο (7,6) από 0 σε 1
      0 0 0 0 0 1 1 0];

% Υπολογισμός των ιδιοτιμών 
eig_values1 = eig(A1);
lambda1_max = max(eig_values1);
lambda1_min = min(eig_values1);

chi_lower1 = 1 + lambda1_max / abs(lambda1_min);
fprintf('Κάτω φράγμα χρωματικού αριθμού για το πρώτο γράφημα: %.2f\n', chi_lower1);

% ---- Επαναλάβετε για το δεύτερο γράφημα ----
A2 = [0 1 1 0 0 0 0 0 0;
      1 0 1 0 0 0 0 1 0;
      1 1 0 0 1 0 0 0 1;
      0 0 0 0 1 1 1 0 0;
      0 0 1 1 0 0 1 0 0;
      0 0 0 1 0 0 1 0 0;
      0 0 0 1 1 1 0 0 1;
      0 1 0 0 0 0 0 0 1;
      0 0 1 0 0 0 0 1 0];

eig_values2 = eig(A2);
lambda2_max = max(eig_values2);
lambda2_min = min(eig_values2);

chi_lower2 = 1 + lambda2_max / abs(lambda2_min);
fprintf('Κάτω φράγμα χρωματικού αριθμού για το δεύτερο γράφημα: %.2f\n', chi_lower2);

% Χρωματισμός του πρώτου γραφήματος
chi_req_1 = ceil(chi_lower1);
colors1 = dsatur_coloring(A1);
fprintf('Χρωματισμός για το πρώτο γράφημα: ');
disp(colors1);
fprintf('Αριθμός χρωμάτων που χρησιμοποιήθηκαν: %d (Αναμενόταν >= %d)\n\n', max(colors1), chi_req_1);

% Χρωματισμός του δεύτερου γραφήματος
chi_req_2 = ceil(chi_lower2);
colors2 = dsatur_coloring(A2);
fprintf('Χρωματισμός για το δεύτερο γράφημα: ');
disp(colors2);
fprintf('Αριθμός χρωμάτων που χρησιμοποιήθηκαν: %d (Αναμενόταν >= %d)\n', max(colors2), chi_req_2);