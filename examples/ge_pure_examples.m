% ================================================================
%   entcalc (MATLAB): ge_pure usage and accuracy discussion
% ================================================================
disp('================================================================');
disp('  entcalc (MATLAB): ge_pure usage and accuracy discussion       ');
disp('================================================================');
disp('This script demonstrates the usage of the ge_pure function.');
disp('It calculates the geometric entanglement for pure states,');
disp('returning both the lower bound and the estimation error caused');
disp('by relaxing the separability condition to the PPT condition.');
disp(' ');
disp('NOTE: The first value is a rigorous lower bound, not just an');
disp('estimate. The true value of the geometric entanglement lies');
disp('exactly at the distance ''error'' from the lower bound (up to');
disp('the numerical accuracy of the underlying SDP solver).');
disp('================================================================');
fprintf('\n');

% ==========================================
% Example 1: 3-Qubit GHZ State
% ==========================================
disp('--- Example 1: 3-Qubit GHZ State ---');
ghz_state = [1; 0; 0; 0; 0; 0; 0; 1];
ghz_state = ghz_state / norm(ghz_state);
dims_3q = [2, 2, 2];

disp('Calculating...');
res_3q = ge_pure(ghz_state, dims_3q);
disp(['Result (Lower bound, Error): [ ', num2str(res_3q), ' ]']);

disp('Notice that the second value in the output (the error) is exactly 0.');
disp('This means the lower bound is strictly equal to the true value.');
disp('This is because the optimization involves a bipartite cut where');
disp('one subsystem is a qubit. For qubit-qubit or qubit-qutrit systems,');
disp('the set of separable states is strictly equal to the set of PPT states.');
fprintf('\n');

% ==========================================
% Example 2: 5-Qubit Random Pure States
% ==========================================
disp('--- Example 2: 5-Qubit Random Pure States ---');
disp('Calculating for 5 completely random 5-qubit states (using QETLAB)...');
dims_5q = [2, 2, 2, 2, 2];

for i = 1:5
    % Generate a random 32-dimensional ket (5 qubits)
    x = RandomStateVector(32);
    res_rand = ge_pure(x, dims_5q);
    disp(['State ', num2str(i), ' -> Result (Lower bound, Error): [ ', num2str(res_rand), ' ]']);
end

disp(' ');
disp('For higher dimensions, there is generally a non-zero estimation error');
disp('caused by the geometric difference between the separable and PPT sets.');
disp('This loop demonstrates that the function can efficiently handle general');
disp('multipartite pure states without relying on any special symmetries.');
fprintf('\n');

disp('================================================================');
disp('Test run completed.');