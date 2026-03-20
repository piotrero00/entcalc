% ================================================================
%   entcalc (MATLAB): upperbip usage and parameter tuning
% ================================================================
disp('================================================================');
disp('  entcalc (MATLAB): upperbip usage and parameter tuning         ');
disp('================================================================');
disp('This script demonstrates the usage of the upperbip function,');
disp('which calculates the upper bound of the geometric entanglement');
disp('for bipartite states (or bipartite cuts) using an iterative');
disp('gradient-descent algorithm.');
disp('Proper parameter usage strongly influences performance and accuracy.');
disp('================================================================');
fprintf('\n');

% ==========================================
% 1. State Preparation (Thermal State)
% ==========================================
disp('--- 1. State Preparation ---');
disp('Constructing a 3-qubit Heisenberg XX spin chain thermal state...');

J = 1.0;

% Define Pauli matrices and Identity matrix for a single qubit
sigmax = [0, 1; 1, 0];
sigmay = [0, -1i; 1i, 0];
sigmaz = [1, 0; 0, -1];
qeye2  = [1, 0; 0, 1];

% Constructing operators for a three-qubit system
sy1 = kron(sigmay, kron(qeye2, qeye2));
sy2 = kron(qeye2, kron(sigmay, qeye2));
sy3 = kron(qeye2, kron(qeye2, sigmay));

sx1 = kron(sigmax, kron(qeye2, qeye2));
sx2 = kron(qeye2, kron(sigmax, qeye2));
sx3 = kron(qeye2, kron(qeye2, sigmax));

% Constant Hamiltonian for the XX model
bet = 1.5;
H = -J/2 * (sy1*sy2 + sy2*sy3 + sy3*sy1 + sx1*sx2 + sx2*sx3 + sx3*sx1);
HH = -H * bet;

% Normalize to get the thermal density matrix rho = expm(-beta*H) / Tr(...)
r = expm(HH);
rho = r / trace(r);

% We evaluate a bipartite cut of the 3-qubit system: 2 qubits vs 1 qubit.
% Hence, the dimensions of the cut are 2^2 = 4 and 2^1 = 2.
dims_bip = [4, 2];
fprintf('State prepared successfully. Bipartite cut dimensions: [4, 2].\n\n');

% ==========================================
% Examples 1 & 2: Random Initialization
% ==========================================
disp('--- Examples 1 & 2: Random Initialization ---');
disp('The algorithm depends on an initial separable decomposition,');
disp('which is generated randomly. For that reason, upperbip called');
disp('with the exact same parameters can give slightly different results.');

% Note: arguments are (rho, dims, dec, {custom_init}, r, iteramax)
res_try1 = upperbip(rho, dims_bip, 0, {0}, 64, 1000);
disp(['iteramax=1000 (first try) : ', num2str(res_try1)]);

res_try2 = upperbip(rho, dims_bip, 0, {0}, 64, 1000);
disp(['iteramax=1000 (second try): ', num2str(res_try2)]);
fprintf('\n');

% ==========================================
% Example 3: Tuning the 'dif' Parameter
% ==========================================
disp('--- Example 3: Tuning the ''dif'' parameter (Threshold) ---');
disp('We can control the accuracy and time of computations using ''dif'' and ''iteramax''.');
disp('In each step, the algorithm increases the fidelity between the input ''rho''');
disp('and the separable decomposition. The algorithm stops when the increase');
disp('in fidelity falls below the ''dif'' threshold.');
disp('Setting ''dif'' too loose (e.g., 1e-5) might cause the algorithm to stop');
disp('too early, before reaching the most accurate (tightest) upper bound.');

% Note: 7th argument is the 'dif' threshold
res_dif_loose = upperbip(rho, dims_bip, 0, {0}, 64, 4500, 1e-5);
disp(['dif=1e-5, iteramax=4500: ', num2str(res_dif_loose)]);
fprintf('\n');

% ==========================================
% Example 4: Tuning the 'iteramax' Parameter
% ==========================================
disp('--- Example 4: Tuning the ''iteramax'' parameter ---');
disp('Let''s tighten the ''dif'' threshold to 1e-9.');
disp('''iteramax'' sets the maximal number of iteration steps. If the number');
disp('of iterations exceeds ''iteramax'', the function stops computations.');
disp('Setting ''iteramax'' too low might cause an inaccurate (loose) result');
disp('because the algorithm hasn''t had enough steps to converge.');

res_dif_tight = upperbip(rho, dims_bip, 0, {0}, 64, 4500, 1e-9);
disp(['dif=1e-9, iteramax=4500: ', num2str(res_dif_tight)]);
fprintf('\n');

disp('================================================================');
disp('Test run completed.');