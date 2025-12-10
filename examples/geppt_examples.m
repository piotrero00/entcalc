
% GHZ state definition
ghz_state = [1; 0; 0; 0; 0; 0; 0; 1]; 
ghz_state = ghz_state / norm(ghz_state);

ghz_proj = ghz_state * ghz_state'; % Density matrix (normalized projection)

w_state = [0; 1; 1; 0; 1; 0; 0; 0];
w_state = w_state / norm(w_state);
w_proj = w_state * w_state'; % Density matrix (normalized projection)

% Mixed state density matrix
rho = 0.7 * ghz_proj + 0.3 * w_proj;

%Default parametrs

disp('Default values: ');
disp(geppt(rho, [2, 2, 2]))

%One can change some parameters
disp('sdpaccuracy=10^(-5): ');
disp(geppt(rho, [2, 2, 2],10^(-5)))