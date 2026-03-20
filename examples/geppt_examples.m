% ================================================================
%   entcalc (MATLAB): geppt usage and parameter tuning
% ================================================================
disp('================================================================');
disp('  entcalc (MATLAB): geppt usage and parameter tuning            ');
disp('================================================================');
disp('This script demonstrates how to use the geppt function in MATLAB.');
disp('It computes a quick lower bound based on the PPT criterion.');
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
disp('geppt is generally the fastest lower bound function, though');
disp('it might not be as tight as ge_mixed_gr or gekppt.');

res_default = geppt(rho, dims);
disp(['Result (Lower bound): ', num2str(res_default)]);
fprintf('\n');

% ==========================================
% Example 2: Tuning Accuracy
% ==========================================
disp('--- Example 2: Tuning accuracy (sdpaccuracy=10^-5) ---');
disp('You can change optional parameters, such as the desired accuracy.');
disp('Decreasing the accuracy can speed up the solver if you only need');
disp('a rough estimate of the PPT lower bound.');

% Note: the 3rd argument is sdpaccuracy
res_acc = geppt(rho, dims, 10^(-5));
disp(['Result (Lower bound): ', num2str(res_acc)]);
fprintf('\n');

disp('================================================================');
disp('Test run completed.');