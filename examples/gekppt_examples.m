% ================================================================
%   entcalc (MATLAB): gekppt usage and parameter tuning
% ================================================================
disp('================================================================');
disp('  entcalc (MATLAB): gekppt usage and k-symmetric extensions     ');
disp('================================================================');
disp('This script demonstrates how different values of ''k'' influence');
disp('the performance of the gekppt function.');
disp('Increasing ''k'' provides a tighter lower bound at the cost of');
disp('significantly higher computational resources.');
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
% Example 1: k=1 Extension
% ==========================================
disp('--- Example 1: k = 1 ---');
disp('We calculate the bound with k=1. This is the fastest, but');
disp('the bound might not be the tightest possible.');
k1 = gekppt(rho, dims, 1);
disp(['Result (k=1): ', num2str(k1)]);
fprintf('\n');

% ==========================================
% Example 2: k=2 Extension
% ==========================================
disp('--- Example 2: k = 2 ---');
disp('We increase k to 2. We expect to obtain a better (tighter)');
disp('lower bound, though it takes more time to compute.');
k2 = gekppt(rho, dims, 2);
disp(['Result (k=2): ', num2str(k2)]);
fprintf('\n');

% ==========================================
% Example 3: k=3 Extension
% ==========================================
disp('--- Example 3: k = 3 ---');
disp('We increase k to 3. This is computationally demanding and');
disp('might take a while. The tightness of the bound increases with k.');
k3 = gekppt(rho, dims, 3);
disp(['Result (k=3): ', num2str(k3)]);
fprintf('\n');

% ==========================================
% Summary of Results
% ==========================================
disp('--- Summary of Lower Bounds ---');
disp('Notice how the lower bound potentially increases (improves)');
disp('as we increase the order of the symmetric extension (k).');
disp('Bounds for [k=1, k=2, k=3]:');
disp([k1, k2, k3]);
fprintf('\n');

disp('================================================================');
disp('Test run completed.');