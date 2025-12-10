function [result] = ge_mixed_gr(state, dim,sdpaccuracy,itera)
% Computes a lower bound on the geometric entanglement of a mixed quantum state 
% using a purity–entanglement complementary relation with higher-dimensional matrices.
%
% Parameters
% ----------
% state : matrix
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
% - This function computes a lower bound using the purity–entanglement complementary relation. 
%   It is one of two functions that use this relation. 
%   This function can be much slower but usually provides a better bound.
%
% Examples
% --------
% GHZ = [1; 0; 0; 0; 0; 0; 0; 1] / sqrt(2);
% W = [0; 1; 1; 0; 1; 0; 0; 0] / sqrt(3);
% GHZ = GHZ*GHZ';
% W = W*W';
% rho = 0.6*GHZ + 0.4*W;
% ge_mixed_gr(rho, [2,2,2])
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
    if ~ishermitian(state)
        error('State must be represented by a positive-semidefinite matrix.');
    end
    if abs(trace(state)-1)>10^(-14)
        error('State is not normalized. Trace of the state should be 1.')
    end
    evals = eig(state);
    if any(real(evals) < -1e-10)
        error('State must be represented by a positive-semidefinite matrix.');
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
    if size(state,1) ~= producta || size(state,2) ~= producta
        error('Dimension of the system is not equal to product of dimensions of subsystem');
    end
   
    
    
    %Main part of the function
    rho = sqrtm(state);
    
  
    ket = rho(:);                   %ket is purification of state
    stp=ket*ket';
   
    pro = (producta^2) ;
    
   
    sta = stp;
    
   
    
    cvx_begin sdp
        
        cvx_solver sdpt3
        cvx_solver_settings( 'gaptol', sdpaccuracy, 'inftol', sdpaccuracy,  'steptol', sdpaccuracy,  'maxit', itera)
        variable X(pro, pro) hermitian semidefinite
    
        X>=0;
        for i = 1:n+1
            partialTranspose_i = PartialTranspose(X, i, [producta, dim]);
            partialTranspose_i >= 0;                                            %PPT relaxation of separability constraint
            
        end
        
    
        
        
        partialTrace_X = PartialTrace(X, 2, [producta, producta]);
        %partialTrace_X == eye(producta);
        
          partialTrace_X == eye(producta);
    
        
        maximize trace(X * sta)
    cvx_end

    if strcmp(cvx_status, 'Inaccurate/Solved')
        warning('Solution did not match with desired precision. Solution may be inaccurate.');
        result = [1 - cvx_optval-sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)),0];
    else
        result = [1 - cvx_optval-sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)),1];        
    end    


end