function result = geppt(rho, dim,sdpaccuracy,itera)
% Computes the lower bound of the geometric entanglement of the state rho 
% using a fidelity relation and PPT relaxation.
%
% Parameters
% ----------
% rho : matrix
%     A 2D array representing a mixed state (density matrix).
% dim : vector of int
%     Dimensions of each subsystem. Each entry should be a positive integer.
% sdpaccuracy : float, optional
%     SDP solver precision. Default is 1e-8.
% itera : int, optional
%     Maximal number of iterations in solver. Default is 200.
%
% Returns
% -------
% ge_lower : array
%     A two-element array:
%       ge_lower(1) – float, lower bound of the geometric entanglement of `state`.
%       ge_lower(2) – int, accuracy flag (1 = accurate, 0 = inaccurate).
%
% Notes
% -----
% - This function computes a lower bound using the relation 
%   E(rho) = 1 - max(F(rho, sigma)), 
%   where the maximization is performed over separable states. 
%   In this implementation, the separability condition is relaxed to PPT, 
%   which provides a lower bound.
%
% Examples
% --------
% GHZ = [1; 0; 0; 0; 0; 0; 0; 1] / sqrt(2);
% W = [0; 1; 1; 0; 1; 0; 0; 0] / sqrt(3);
% GHZ = GHZ*GHZ';
% W = W*W';
% rho = 0.6*GHZ + 0.4*W;
% geppt(rho, [2,2,2])
% ans =
%    [0.3025,1.0000]
    if nargin < 3 || isempty(sdpaccuracy)
        sdpaccuracy = 1e-8;
    end
    if nargin < 4 || isempty(itera)
        itera = 200;
    end
    n = length(dim);
    if min(dim) <= 1
        error('Dimension of subsystem must be an integer greater than 1');
    end
    if any(floor(dim) ~= dim)
        error('Dimension of subsystem must be an integer greater than 1');
    end
    if n <= 1
        error('System must be at least bipartite');
    end
    if ~ishermitian(rho)
        error('State must be represented by a positive-semidefinite matrix.');
    end

    
    evals = eig(rho);
    if any(real(evals) < -1e-10)
        error('State must be represented by a positive-semidefinite matrix.');
    end

    if abs(trace(rho)-1)>10^(-14)
        error('State is not normalized. Trace of the state should be 1.')
    end
    if ~isa(sdpaccuracy,'double')
        error('sdpaccuracy must be a float greater than 0');
    end
    if sdpaccuracy >=1 
        error('sdpaccuracy must be a float smaller than 1');
    end
    if sdpaccuracy <= 0
        error('sdpaccuracy must be a float greater than 0');
    end

    
    if floor(itera) ~= itera
        error('itera must be an integer greater than 0');
    end
    if itera <= 0
        error('itera must be an integer greater than 0');
    end

    
    producta = prod(dim);
    if size(rho,1) ~= producta || size(rho,2) ~= producta
        error('Dimension of the system is not equal to product of dimensions of subsystem');
    end
   
    
    
    %Main part of the function
    
    
    pro = prod(dim); 
    
    zer = complex(zeros(pro, pro));       
    I_pro = eye(pro);                    
    
    
    A = [zer, I_pro/2;
         I_pro/2, zer];
    AA = [I_pro, zer;
         zer, zer];
    cvx_begin sdp
        cvx_solver sdpt3
        cvx_solver_settings( 'gaptol', sdpaccuracy, 'inftol', sdpaccuracy,  'steptol', sdpaccuracy,  'maxit', itera)

        variable siigma(pro, pro) hermitian semidefinite
        variable X(pro, pro) complex
        dual variable Y
        dual variable lb
        
        lb : trace(siigma) == 1;
        siigma >= 0;
        
        for i=1:1:n
            sigma_pt = PartialTranspose(siigma, i, dim); 
            sigma_pt >= 0;
        end
       
        B = [rho, X;
             X', siigma];
        
        Y : B >= 0;
        
       
        maximize trace(A*B)
    cvx_end
   % trace(A*B)
   % Y11=Y(1:pro,1:pro);
   % trace(Y11*rho)-lb
   % cvx_optval
   % cvx_optval+sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy))
    if strcmp(cvx_status, 'Inaccurate/Solved')
        warning('Solution did not match with desired precision. Solution may be inaccurate.');
        
        result = [1 - (cvx_optval+sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)))^2,0];
    else
        result = [1 - (cvx_optval+sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)))^2,1];
    end
end