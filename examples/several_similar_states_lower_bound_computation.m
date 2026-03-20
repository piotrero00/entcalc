% ================================================================
%   entcalc (MATLAB): uppermul and 'warm start' optimization
% ================================================================
disp('================================================================');
disp('  entcalc (MATLAB): uppermul and ''warm start'' optimization    ');
disp('================================================================');
disp('This script demonstrates an advanced usage of the uppermul');
disp('function. It compares a naive computation approach with an');
disp('accelerated ''warm start'' approach for similar states.');
disp('By passing the decomposition components (probabilities and kets)');
disp('from a previous state, we can drastically reduce iterations.');
disp(' ');
disp('[WARNING] Execution of this script might take a significant');
disp('amount of time. Please be patient.');
disp('================================================================');
fprintf('\n');

% ==========================================
% 1. Operator Setup (Global)
% ==========================================
disp('--- 1. Initializing Hamiltonian and Thermal States ---');
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
H = -J/2 * (sy1*sy2 + sy2*sy3 + sy3*sy1 + sx1*sx2 + sx2*sx3 + sx3*sx1);

disp('Hamiltonian constructed. Starting computations...');
fprintf('\n');

% ==========================================
% 2. Naive Approach: Random Initialization
% ==========================================
disp('--- 2. Starting Naive Computations (Random Start) ---');
disp('Computing upper bounds for varying temperatures (beta from 0.60 to 0.95).');
disp('Starting from a random decomposition requires a high iteramax (4000).');

results_naive = [];
tic; % Start timer for naive approach

for i = 60:5:95 
    bet = i/100;
    
    HH = -H * bet;
    r = expm(HH);
    state = r / trace(r);
        
    % 4th argument {0} means no custom initialization
    res = uppermul(state, [2, 2, 2], 0, {0}, 64, 4000);
    
    results_naive = [results_naive, res]; % Append the upper bound
    disp(['Beta: ', num2str(bet, '%.2f'), ' | Upper bound: ', num2str(res, '%.9f')]);
end

time_naive = toc; % Stop timer
disp(['-> Naive approach took: ', num2str(time_naive, '%.2f'), ' seconds.']);
fprintf('\n');

% ==========================================
% 3. Accelerated Approach: Warm Start
% ==========================================
disp('--- 3. Starting Accelerated Computations (Warm Start) ---');
disp('We use the decomposition of the previous state as the initial seed');
disp('for the next state. This allows us to slash iteramax from 4000 to 500!');

results_acc = [];
tic; % Start timer for accelerated approach

% Initial reference state (bet=0.6)
bet_init = 0.6;
HH = -H * bet_init;
r = expm(HH);
state = r / trace(r);

disp(['Computing initial seed for Beta: ', num2str(bet_init, '%.2f'), '...']);
% First computation must be done from random decomposition.
% We set dec=1 (3rd arg) to get the decomposition for the next step.
res_cell = uppermul(state, [2, 2, 2], 1, {0}, 64, 4000);

% res_cell{1} is the upper bound, res_cell{2} are probs, res_cell{3} are kets
results_acc = [results_acc, res_cell{1}]; 
disp(['Beta: ', num2str(bet_init, '%.2f'), ' | Upper bound (seed) : ', num2str(res_cell{1}, '%.9f')]);

% Iterative loop using previous results as seeds
for i = 65:5:95 
    bet = i/100;
    
    % Recalculate state for new beta
    HH = -H * bet;
    r = expm(HH);
    state = r / trace(r);
        
    % Now we start with the decomposition from the previous step.
    % We pass the previous decomposition (res_cell{2}, res_cell{3}) via the 4th parameter.
    % Format: {1 (flag for custom init), probabilities, kets}
    res_cell = uppermul(state, [2, 2, 2], 1, {1, res_cell{2}, res_cell{3}}, 64, 500);

    results_acc = [results_acc, res_cell{1}]; % Append the new upper bound
    disp(['Beta: ', num2str(bet, '%.2f'), ' | Upper bound (accel): ', num2str(res_cell{1}, '%.9f')]);
end

time_acc = toc; % Stop timer
disp(['-> Accelerated approach took: ', num2str(time_acc, '%.2f'), ' seconds.']);

fprintf('\n');
disp('================================================================');
disp('Test run completed. Notice the massive difference in execution time!');