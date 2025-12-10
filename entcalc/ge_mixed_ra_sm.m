function [result] = ge_mixed_ra_sm(kets,ps, dim,sdpaccuracy,itera)
% Computes a lower bound on the geometric entanglement of a mixed quantum state 
% using a purity–entanglement complementary relation with higher-dimensional matrices.
%Contrary to ge_mixed_sm this function accepts input as a ensemble of orthogonal kets not a density matrix.
%Because of this it can find a more efficient purification.
%
% Parameters
% ----------
%kets : cell
%    A cell of kets (arrays representing a pure state). Kets must be orthogonal.
%ps : vector of floats.
%    A vector of non-negative floats representing probability of given ket in state ensemble. 
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
% ge_mixed_ra_sm({GHZ,W},[0.6,0.4], [2,2,2])
% ans =
%    [0.3025,1.0000]
    if nargin < 4 || isempty(sdpaccuracy)
        sdpaccuracy = 1e-8;
    end
    if nargin < 5 || isempty(itera)
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
    
    rank=length(kets);
    for i = 1:rank
        for j = i+1:rank
            overlap = kets{i}' * kets{j};   % <i|j>
            if abs(overlap) > 1e-14
                error('Kets should be orthogonal');
            end
        end
    end
    

    nk = length(kets{1});  

    for i = 1:rank
        
        if abs(norm(kets{i}) - 1) >= 1e-14
            error('States should be normalized');
        end
        
        
        if length(kets{i}) ~= nk
            error('States should have the same dimensions');
        end
    end
    producta = prod(dim);
    if nk ~= producta
        error('Dimension of the system is not equal to product of dimensions of subsystem');
    end
    
    if abs(sum(ps) - 1) > 1e-14
        error('Probabilities do not sum up to 1');
    end
    
    if min(ps) <= -1e-14
        error('Probabilities should be greater than 0');
    end
    
    
    
   %Main part of the function
    
    
    stb = zeros(producta*rank, 1);  
    for i = 1:rank                   %We construct purification
        e_i = basis(rank, i);        
        term = sqrt(ps(i)) * kron(e_i, kets{i});  
        stb = stb + term;
    end

    
                 %stb is purification of state
    stc=stb*stb';
    

    sta=PartialTrace(stc,n+1,[rank,dim]);
    dim(end)=[];
    product=prod(dim);
    pro = rank*product;
    
   
    
    
   
    
    cvx_begin sdp
        
        cvx_solver sdpt3
        cvx_solver_settings( 'gaptol', sdpaccuracy, 'inftol', sdpaccuracy,  'steptol', sdpaccuracy,  'maxit', itera)
        variable X(pro, pro) hermitian semidefinite
    
        X>=0;
        for i = 1:n
            partialTranspose_i = PartialTranspose(X, i, [rank, dim]);
            partialTranspose_i >= 0;                                            %PPT relaxation of separability constraint
            
        end
        
    
        
        
        partialTrace_X = PartialTrace(X, 2, [rank, product]);
        
        dual variable Z
        Z : partialTrace_X == eye(rank);
    
        
        maximize trace(X * sta)
    cvx_end

    if strcmp(cvx_status, 'Inaccurate/Solved')
        warning('Solution did not match with desired precision. Solution may be inaccurate.');
        result = [1 - cvx_optval-sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)),0];
    else
        result = [1 - cvx_optval-sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)),1];        
    end    


end