% ================================================================
%   entcalc (MATLAB): uppermul usage and accuracy discussion
% ================================================================
disp('================================================================');
disp('  entcalc (MATLAB): uppermul usage and accuracy discussion      ');
disp('================================================================');
disp('This script demonstrates the usage of the uppermul function');
disp('with optional parameters. We discuss the numerical accuracy');
disp('of computations and how to extract the separable decomposition.');
disp(' ');
disp('[WARNING] Testing parameters on a fully random state requires');
disp('more iterations to converge. This script might take a moment.');
disp('================================================================');
fprintf('\n');

% ==========================================
% 1. State Preparation (GHZ Pure State)
% ==========================================
% GHZ state definition
ghz_state = [1; 0; 0; 0; 0; 0; 0; 1]; 
ghz_state = ghz_state / norm(ghz_state);
ghz_proj = ghz_state * ghz_state'; % Density matrix (normalized projection)

% ==========================================
% Example 1: Numerical Precision (Machine Epsilon)
% ==========================================
disp('--- Example 1: Numerical Precision and Machine Epsilon ---');
disp('We run uppermul 20 times on the exact same GHZ state.');
disp('Calculating...');

results = [];
for i = 1:20
    % 4th argument {0} implies no custom warm-start initialization
    res = uppermul(ghz_proj, [2, 2, 2], 0, {0});
    results = [results, res];
end

disp(['Minimum result: ', num2str(min(results), '%.10f')]);
disp(['Maximum result: ', num2str(max(results), '%.10f')]);
fprintf('\n');

disp('It is known that the geometric entanglement of the GHZ state is 0.5.');
disp('We obtain results that might be around 10^-9 greater OR lower.');
disp('Results slightly lower than 0.5 might be surprising, because these');
disp('functions should output an UPPER bound. This happens due to machine');
disp('precision. Since we are taking square roots of machine precision (10^-16),');
disp('we might encounter a ~10^-8 error in the results. One needs to keep');
disp('in mind that the strict upper bound can be this amount higher than the output.');
fprintf('\n');

% ==========================================
% Example 2: Extracting Separable Decomposition
% ==========================================
disp('--- Example 2: Verifying the Decomposition (dec=1) ---');
disp('If you have doubts about the result, you can check the correctness');
disp('of the separable decomposition using the optional argument dec=1.');
disp('The output becomes a cell array containing the bound, probabilities, and kets.');

res_dec = uppermul(ghz_proj, [2, 2, 2], 1, {0});
disp(['Upper bound: ', num2str(res_dec{1})]);
disp('Probability distribution and kets are stored in res_dec{2} and res_dec{3}.');
fprintf('\n');

% ==========================================
% 2. State Preparation (GHZ-W Mixture)
% ==========================================
% W state definition
w_state = [0; 1; 1; 0; 1; 0; 0; 0];
w_state = w_state / norm(w_state);
w_proj = w_state * w_state'; % Density matrix 

% Mixed density matrix state
rho = 0.6 * ghz_proj + 0.4 * w_proj;

disp('================================================================');
disp('  Parameter Tuning for Complex States                           ');
disp('================================================================');

disp('--- Example 3: Default values (GHZ-W mixture) ---');
b = uppermul(rho, [2, 2, 2], 0, {0});
disp(['Default values result for simple rank-2 state: ', num2str(b)]);
fprintf('\n');

% ==========================================
% 3. Complex State for Parameter Tuning
% ==========================================
disp('--- Testing parameters on a complex state ---');
disp('The GHZ-W mixture above is a simple, rank-2 state. The algorithm finds');
disp('the global optimum so easily that parameter tweaking doesn''t change much.');
disp('To truly see the impact of parameters, we generate a random 3-qubit state.');

% Generate a completely random, full-rank density matrix using QETLAB
rho_hard = RandomDensityMatrix(8);

% ==========================================
% Example 4: The 'r' Parameter
% ==========================================
disp('--- Example 4: The ''r'' Parameter (r=10) ---');
disp('By Caratheodory''s theorem, the optimal separable decomposition has r <= d^2.');
disp('We can accelerate computations by taking a smaller r (e.g., r=10). Keep in mind');
disp('that we might decrease the accuracy by doing this.');

res_r4 = uppermul(rho_hard, [2, 2, 2], 0, {0}, 10);
disp(['Result (r=10): ', num2str(res_r4)]);
fprintf('\n');

% ==========================================
% Example 5: The 'sepitera' Parameter
% ==========================================
disp('--- Example 5: The ''sepitera'' Parameter (sepitera=25) ---');
disp('By increasing ''sepitera'', we can increase the number of outer iterations');
disp('before convergence. But be careful: setting ''sepitera'' too high slows');
disp('down computations without significantly increasing accuracy.');

% args: state, dims, dec, init, r, iteramax, dif, sepitera
res_sep = uppermul(rho_hard, [2, 2, 2], 0, {0}, 64, 1000, 10^(-9), 25);
disp(['Result (sepitera=25): ', num2str(res_sep)]);
fprintf('\n');

% ==========================================
% Example 6: The 'dif' Parameter
% ==========================================
disp('--- Example 6: The ''dif'' Parameter (dif=10^-5) ---');
disp('A larger ''dif'' makes computations faster, but the result might not be');
disp('accurate. The function terminates if the fidelity update between subsequent');
disp('iterations is smaller than ''dif''.');

res_dif = uppermul(rho_hard, [2, 2, 2], 0, {0}, 64, 1000, 10^(-5));
disp(['Result (dif=10^-5): ', num2str(res_dif)]);
fprintf('\n');

% ==========================================
% Example 7: The 'iteramax' Parameter
% ==========================================
disp('--- Example 7: The ''iteramax'' Parameter (10 vs 1000) ---');
disp('''iteramax'' determines the maximal number of iterations. It prevents');
disp('too long computations. Setting it too low results in an imprecise bound.');

res_it_low = uppermul(rho_hard, [2, 2, 2], 0, {0}, 64, 10);
disp(['Result (iteramax=10 - very fast, inaccurate): ', num2str(res_it_low)]);

res_it_high = uppermul(rho_hard, [2, 2, 2], 0, {0}, 64, 1000);
disp(['Result (iteramax=1000 - slower, highly accurate): ', num2str(res_it_high)]);
fprintf('\n');

disp('================================================================');
disp('Test run completed.');