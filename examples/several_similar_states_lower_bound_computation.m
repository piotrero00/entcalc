%
% We demonstrate how to use uppersame to accelerate computations of
% entanglement upper bound for similar states. Assume we want to compute
% the upper bound for some spin chains with changing inverse temperature.
%

% Define set-up

J = 1.0;

% Define Pauli matrices and Identity matrix for a single qubit (dimension 2)
sigmax = [0, 1; 1, 0];
sigmay = [0, -1i; 1i, 0];
sigmaz = [1, 0; 0, -1];
qeye2 = [1, 0; 0, 1];

% Constructing operators for a three-qubit system (dimension 8)
% Equivalent to qutip.tensor operations

sy1 = kron(sigmay, kron(qeye2, qeye2));
sy2 = kron(qeye2, kron(sigmay, qeye2));
sy3 = kron(qeye2, kron(qeye2, sigmay));

sx1 = kron(sigmax, kron(qeye2, qeye2));
sx2 = kron(qeye2, kron(sigmax, qeye2));
sx3 = kron(qeye2, kron(qeye2, sigmax));

sz1 = kron(sigmaz, kron(qeye2, qeye2));
sz2 = kron(qeye2, kron(sigmaz, qeye2));
sz3 = kron(qeye2, kron(qeye2, sigmaz));

% Initial inverse temperature
bet = 0.6;
% Hamiltonian
H = -J/2 * (sy1*sy2 + sy2*sy3 + sy3*sy1 + sx1*sx2 + sx2*sx3 + sx3*sx1);
HH = -H * bet;

% Compute and normalize initial state (density matrix)
r = expm(HH);
state = r / trace(r);

% Initialize results storage
results = [];

% Perform the first calculation to get an initial decomposition (not strictly necessary here, 
% as it's overwritten by the loop, but kept for structural consistency)
res = uppermul(state, [2, 2, 2],0,{0},64,4000);

% --- Loop 1: Standard Calculation (Starting from random decomposition) ---
for i = 60:5:95 
    bet = i/100;
    
    
    H = -J/2 * (sy1*sy2 + sy2*sy3 + sy3*sy1 + sx1*sx2 + sx2*sx3 + sx3*sx1);
    HH = -H * bet;
    r = expm(HH);
    state = r / trace(r);
        
    
    res = uppermul(state, [2, 2, 2],0,{0},64,4000);
    
    results = [results, res]; % Append the upper bound
end

%
% We can accelerate computations by starting from state that was output of the previous state.
% Doing that we do not have random decomposition at the beginning and we converge faster.
% We can put lower iteration number.
%

% --- Loop 2: Accelerated Calculation (Using previous decomposition) ---

% Reset initial state (bet=0.6)
bet = 0.6;
H = -J/2 * (sy1*sy2 + sy2*sy3 + sy3*sy1 + sx1*sx2 + sx2*sx3 + sx3*sx1);
HH = -H * bet;
r = expm(HH);
state = r / trace(r);

result2 = [];

% First computation must be done from random decomposition.
% We need 'dec=1' to get the decomposition for the next step.
res = uppermul(state, [2, 2, 2],1,{0},64,4000);

result2 = [result2, res{1}]; % First output is the upper bound (res{1} in MATLAB cell array)

% The decomposition is stored in res{2} (probabilities) and res{3} (kets/states)

for i = 65:5:95 
    bet = i/100;
    
    % Recalculate state for new beta
    H = -J/2 * (sy1*sy2 + sy2*sy3 + sy3*sy1 + sx1*sx2 + sx2*sx3 + sx3*sx1);
    HH = -H * bet;
    r = expm(HH);
    state = r / trace(r);
        
    %
    % Now we will start with decomposition, which optimized previous computations.
    % Since change in state is small, we expect the output to converge faster.
    %
    
    % Secend output is probability distribution, third output is state ensamble
    % We pass the previous decomposition (res{2}, res{3}) via 'qs' and 'sqs' parameters.
    res = uppermul(state, [2, 2, 2],1,{1,res{2},res{3}},64,500);

    result2 = [result2, res{1}]; % Append the new upper bound

    % The next iteration will use the decomposition stored in 'res{2} and res{3}' from this step.
end


%
% The output res{3} consists of a matrix whose columns are the subsequent kets of the separable decomposition.
% res{2} outputs probabilities corresponding to these kets.
%