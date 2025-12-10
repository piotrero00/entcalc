% MATLAB does not require explicit imports like Python, but assumes necessary
% functions (e.g., upperbip_matlab) are available in the path.

J = 1.0;

% Define Pauli matrices and Identity matrix for a single qubit (dimension 2)
% Equivalent to qutip.sigmax/y/z and qutip.qeye(2)
sigmax = [0, 1; 1, 0];
sigmay = [0, -1i; 1i, 0];
sigmaz = [1, 0; 0, -1];
qeye2 = [1, 0; 0, 1];

% Constructing operators for a three-qubit system (dimension 8)
% Use Kronecker product (kron) for tensor products.

% Sigma Y operators (SY1, SY2, SY3)
sy1 = kron(sigmay, kron(qeye2, qeye2));
sy2 = kron(qeye2, kron(sigmay, qeye2));
sy3 = kron(qeye2, kron(qeye2, sigmay));

% Sigma X operators (SX1, SX2, SX3)
sx1 = kron(sigmax, kron(qeye2, qeye2));
sx2 = kron(qeye2, kron(sigmax, qeye2));
sx3 = kron(qeye2, kron(qeye2, sigmax));

% Sigma Z operators (SZ1, SZ2, SZ3)
sz1 = kron(sigmaz, kron(qeye2, qeye2));
sz2 = kron(qeye2, kron(sigmaz, qeye2));
sz3 = kron(qeye2, kron(qeye2, sigmaz));


% Define inverse temperature parameter and Hamiltonian
bet = 1.5;

% Define the Hamiltonian H for the three-qubit system (J/2 * sum of couplings)
H = -J/2 * (sy1*sy2 + sy2*sy3 + sy3*sy1 + sx1*sx2 + sx2*sx3 + sx3*sx1);

% Compute -H*bet (for exponentiation)
HH = -H * bet;

% Compute the unnormalized density matrix r = exp(-beta*H)
% Equivalent to r = HH.expm()
r = expm(HH);

% Normalize to get the thermal density matrix rho = r / Tr(r)
rho = r / trace(r);


%
% We demonstrate different usage of parameters for upperbip.
% Proper parameter usage influences the performance of the functions.
%


disp('iteramax=1000 first try: ');
 disp(upperbip(rho, [4, 2],0,{0},64 ,1000));

%
% The algorithm depends on initial separable decomposition, which is random.
% For that reason, upperbip with the same parameters can give different result.
%
disp('iteramax=1000 second try: ');
disp(upperbip(rho, [4, 2],0,{0},64 ,1000));

%
% We can control accuracy and time of computations, by iteramax, dif.


disp('dif=10^(-5),iteramax=4500: ');
disp(upperbip(rho, [4, 2],0,{0},64,4500,1e-5))
%
% Too little dif might cause the algorithm to stop before reaching an accurate bound.
% In each algorithm step we increase the fidelity of input rho and separable decomposition.
% The algorithm stops when the increase in fidelity falls below dif.
% For that reason, setting too little dif might cause too early stop.
%

disp('dif=10**(-9),iteramax=4500: ');
disp(upperbip(rho, [4, 2],0,{0},64,4500,1e-9))
%
% iteramax sets the maximal number of iteration steps. If the number of iterations exceeds iteramax,
% the function stops computations. Setting too little iteramax might cause an inaccurate result.
%