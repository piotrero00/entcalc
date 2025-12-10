
ghz_state = [1; 0; 0; 0; 0; 0; 0; 1]; % Example of an 8-dimensional vector
ghz_state = ghz_state / norm(ghz_state); % Normalization
ghz_proj = ghz_state * ghz_state'; % Density matrix (projection)

% Section Title: Function Usage and Accuracy Discussion
%
% In this script we show some examples of usage of uppersame with optional parameters.
% We also discuss accuracy of computations.
%

results = [];

for i = 1:20
    
    res = uppermul(ghz_proj, [2, 2, 2],0,{0});
    results = [results, res]; % Appending the result to the results vector
end

disp([min(results), max(results)]) % Displaying the minimum and maximum results

%
% It is known that the geometric entanglement of GHZ state is 0.5. We obtain results
% around 10**(-9) greater and lower.
% Results slightly lower than 0.5 might be surprising, because uppersame should
% output the upper bound. It happens because of machine precision accuracy.
% Since we are taking square roots of machine precision 10**(-16), we might encounter
% 10**(-8) error in results. One needs to keep in mind that accuracy is around 10**(-8)
% and the upper bound can be this amount higher than the output of the uppersame.
% If one has some doubts about the result, one can check the correctness of the separable
% decomposition. If one uses optional argument dec=1, the separable decomposition is
% also returned. One can then check that it is a true separable decomposition.
%

% Using the optional 'dec' argument

res_dec = uppermul(ghz_proj, [2, 2, 2], 1, {0});

disp('Upper bound: ')
disp(res_dec{1})

% Another state: The W state

w_state = [0; 1; 1; 0; 1; 0; 0; 0]; % Example vector
w_state = w_state / norm(w_state); % Normalization
w_proj = w_state * w_state'; % Density matrix

% Mixed density matrix state
rho = 0.6 * ghz_proj + 0.4 * w_proj;
disp("Results for GHZ-W mixture:")
disp("Default values result")
b=uppermul(rho, [2, 2, 2],0,{0});
disp(b)

% Using the optional 'r' argument
disp("r=4")
res_r10 = uppermul(rho, [2, 2, 2], 0, {0},4);
disp(res_r10)

%
% We accelerate computations by taking smaller r. Keep in mind that we can decrease
% accuracy by doing this.
% By Carathéodory's theorem we are guaranteed that optimal separable decomposition
% has r=d**(2), where d is the number of dimensions of rho.
% Parameters 'sepitera', 'dif', 'iteramax' determine accuracy and speed of computations.
%

% Using the optional 'sepitera' argument

res_sepitera = uppermul(rho, [2, 2, 2], 0,{0}, 64, 1000,10^(-9),25);
disp("sepitera=25")
disp(res_sepitera)

%
% Increasing sepitera we can reduce the number of iterations before convergence.
% But be careful, too great sepitera slows down computations without increasing accuracy.
%

% Using the optional 'dif' argument

res_dif = uppermul(rho, [2, 2, 2], 0,{0}, 64, 10,10^(-5));
disp('dif=10**(-5): ')
disp(res_dif)

%
% Smaller dif makes computations faster, but the result might not be accurate.
% The function terminates if the fidelity update between subsequent iterations is smaller than dif.
% If the function converges slowly, putting too little dif might result in an unprecise bound.
%

% Using the optional 'iteramax' argument

res_iteramax = uppermul(rho, [2, 2, 2], 0,{0}, 64, 10);
disp("iteramax=10")
disp(res_iteramax)

% Note that with dec=0 (default) the output is a float number, with dec=1 it is a list (cell array in MATLAB) with 3 elements.

%
% iteramax determines the maximal number of iterations. It prevents too long computations.
% Together with dif it can be used to manipulate the precision of the computations.
%