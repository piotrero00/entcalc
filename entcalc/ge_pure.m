function sol = ge_pure(rho, dim)
%Returns the geometric entanglement of a pure state `rho`.
%Returns the geometric entanglement of a pure state `rho`.
%
% `rho` should be given in ket form. `dim` is a vector of dimensions
% of each subsystem. `n` is the number of subsystems.
%
% Parameters
% ----------
% rho : vector
%     A 1D array representing a pure state in ket form.
% dim : vector of int
%     Dimensions of each subsystem. Each entry should be a positive integer.
%
% Returns
% -------
% cell
%     A cell containing:
%       - The geometric entanglement of `rho` (float).
%       - The numerical error of the result (float).
%
% Notes
% -----
% - For qubit-qubit and qubit-qutrit systems, this function gives exact solutions.
%   The only error is the numerical error connected with SDP precision.
%
% Examples
% --------
% GHZ = [1; 0; 0; 0; 0; 0; 0; 1] / sqrt(2);
% ge_pure(GHZ, [2,2,2])
% ans =
%     {0.500000, 0.0}   % The only error is connected with SDP precision

    n=length(dim);
                           %Checking correctness of input
    if min(dim) <= 1
        error('Dimension of subsystem must be an integer greater than 1');
    end
    if n < 3
        error('State must be at least tripartite');
    end
    if any(floor(dim) ~= dim)
        error('Dimension of subsystem must be an integer greater than 1');
    end
    
    
    producta = prod(dim);
    if length(rho) ~= producta
        error('Dimension of the system is not equal to product of dimensions of subsystem');
    end
    
    
    if abs(1 - norm(rho)) > 1e-8
        warning('Warning! State is not normalized. The result in that scenario is not geometric entanglement.');
    end
    
    
    sz = size(rho);
    if ~(isvector(rho) && (sz(2) == 1))  
        error('State should be given in ket form.');
    end
    [~, max_index] = max(dim);

    
    li = 1:n;
    li(max_index) = [];
    dim_reduced = dim;                  %We trace out the largest subsystemphi
    dim_reduced(max_index) = [];
    
   
    traceout = setdiff(1:n, li);

    
    
    rh = PartialTrace(rho*rho', traceout, dim);
   
    product = prod(dim_reduced);

   
    cvx_begin sdp quiet
        cvx_solver sdpt3
        cvx_precision high
        variable X(product,product) hermitian semidefinite

        maximize( real(trace(rh * X)) )
        subject to
            trace(X) == 1;
            X>=0

           
            for i = 1:length(dim_reduced)
                PartialTranspose(X, i,dim_reduced) >= 0;
            end
    cvx_end

   
    ww = real(eig(X));
    maxww = max(ww);
    ep = abs(maxww - 1);

    if abs(maxww - 1) > 1e-8 && product > 6
        warning('Warning! entcalc could not find good approximation bound.');
    end

   
    if product <= 6 && n == 3
        accura = 0;
    else
        if n == 3
            accura = 4*sqrt(ep);
        elseif n == 4
            accura = 8*sqrt(ep);
        else
            accura = 69; 
        end
    end
    if accura > 0.01
        if abs(sum(X(:)) - trace(X)) < 1e-10
            sig_dia = diag(diag(X));       %diagonal states are certainly separable.
            acc = trace(sig_dia * rh);
            accura = abs(cvx_optval - acc);
        end
    end
    global cvx_solver_output
    sol = [1 - cvx_optval,accura, cvx_solver_output];
end
