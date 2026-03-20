% ================================================================
%   entcalc (MATLAB): ge_mixed_gr usage and parameter tuning
% ================================================================
disp('================================================================');
disp('  entcalc (MATLAB): ge_mixed_gr usage and parameter tuning      ');
disp('================================================================');
disp('This script demonstrates how to use the ge_mixed_gr function');
disp('in MATLAB to compute the lower bound of geometric entanglement.');
disp('================================================================');
fprintf('\n');

% ==========================================
% 1. State Preparation
% ==========================================
disp('--- 1. State Preparation ---');
disp('Creating a 3-qubit mixture: 70% GHZ and 30% W state...');

% GHZ state definition
ghz_state = [1; 0; 0; 0; 0; 0; 0; 1]; 
ghz_state = ghz_state / norm(ghz_state);
ghz_proj = ghz_state * ghz_state'; % Density matrix (normalized projection)

% W state definition
w_state = [0; 1; 1; 0; 1; 0; 0; 0];
w_state = w_state / norm(w_state);
w_proj = w_state * w_state'; % Density matrix (normalized projection)

% Mixed state density matrix
rho = 0.7 * ghz_proj + 0.3 * w_proj;

% Subsystem dimensions for 3 qubits
dims = [2, 2, 2];
fprintf('State prepared successfully.\n\n');

% ==========================================
% Example 1: Default Parameters
% ==========================================
disp('--- Example 1: Default values ---');
disp('Calculating with default parameters...');
disp('We usually find an accurate solution with these settings.');

res_default = ge_mixed_gr(rho, dims);
disp(['Result (Lower bound): ', num2str(res_default)]);
fprintf('\n');

% ==========================================
% Example 2: Decreasing Required Accuracy
% ==========================================
disp('--- Example 2: Decreasing accuracy (sdpaccuracy=10^-5) ---');
disp('One can change some parameters, such as desired accuracy.');
disp('The true lower bound lies within the sdpaccuracy range from the result.');
disp('If you need to quickly estimate the bound, decreasing the required');
disp('accuracy is a good option.');

res_acc = ge_mixed_gr(rho, dims, 10^(-5));
disp(['Result (Lower bound): ', num2str(res_acc)]);
fprintf('\n');

disp('================================================================');
disp('Test run completed.');