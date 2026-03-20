% ================================================================
%   entcalc (MATLAB): ge_mixed_ra_gr usage and parameter tuning
% ================================================================
disp('================================================================');
disp('  entcalc (MATLAB): ge_mixed_ra_gr usage and parameter tuning   ');
disp('================================================================');
disp('This script demonstrates how to use the ge_mixed_ra_gr function');
disp('in MATLAB. Unlike ge_mixed_gr, this function takes a cell array');
disp('of pure state vectors and their corresponding probabilities.');
disp('Thanks to this, it computes the lower bound more efficiently.');
disp('================================================================');
fprintf('\n');

% ==========================================
% 1. State Preparation (Ensemble)
% ==========================================
disp('--- 1. State Preparation (Ensemble) ---');
disp('Creating an ensemble of 3-qubit states: 70% GHZ and 30% W state...');

% GHZ state definition (column vector)
ghz_state = [1; 0; 0; 0; 0; 0; 0; 1]; 
ghz_state = ghz_state / norm(ghz_state);

% W state definition (column vector)
w_state = [0; 1; 1; 0; 1; 0; 0; 0];
w_state = w_state / norm(w_state);

% Define the ensemble: cell array of states and array of probabilities
states = {ghz_state, w_state};
probs = [0.7, 0.3];

% Subsystem dimensions for 3 qubits
dims = [2, 2, 2];
fprintf('Ensemble prepared successfully.\n\n');

% ==========================================
% Example 1: Default Parameters
% ==========================================
disp('--- Example 1: Default values ---');
disp('Calculating with default parameters...');
disp('We usually find an accurate solution with these settings.');

res_default = ge_mixed_ra_gr(states, probs, dims);
disp(['Result (Lower bound): ', num2str(res_default)]);
fprintf('\n');

% ==========================================
% Example 2: Decreasing Required Accuracy
% ==========================================
disp('--- Example 2: Decreasing accuracy (sdpaccuracy=10^-5) ---');
disp('One can change optional parameters, such as desired accuracy.');
disp('The true lower bound lies within the sdpaccuracy range from the result.');
disp('If you need to quickly estimate the bound, decreasing the required');
disp('accuracy is a good option.');

% Note: the 4th argument is sdpaccuracy
res_acc = ge_mixed_ra_gr(states, probs, dims, 10^(-5));
disp(['Result (Lower bound): ', num2str(res_acc)]);
fprintf('\n');

disp('================================================================');
disp('Test run completed.');