function result = gekppt(rho, dim, k,sdpaccuracy,itera)
% Computes the lower bound of the geometric entanglement of the state `rho` 
% using fidelity relation, PPT relaxation, and k-symmetric extensions.
%
% Parameters
% ----------
% rho : matrix
%     A 2D array representing a mixed state (density matrix).
% dim : vector of int
%     Dimensions of each subsystem. Each entry should be a positive integer.
% k : integer
%     Number of k-symmetric extensions.
%     For a bipartite state ρ^(AB), this specifies the number of copies of subsystem B.
%     For example, if k = 2, the extended state is ρ^(ABB').
%     For multipartite states, k defines the number of additional copies 
%     of each subsystem beyond the first one.
%     For example, for a tripartite state ρ^(ABC) and k = 4,
%     the extended state is ρ^(ABCB'C'B''C'').
% sdpaccuracy : float, optional
%     SDP precision. Default is 1e-8.
% itera : int, optional
%     Maximal number of iterations in solver. Default is 200.
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
%   In this implementation, the separability condition is relaxed 
%   to PPT together with k-symmetric extendibility, 
%   which provides a lower bound.
%
% Examples
% --------
% GHZ = [1; 0; 0; 0; 0; 0; 0; 1] / sqrt(2);
% W = [0; 1; 1; 0; 1; 0; 0; 0] / sqrt(3);
% GHZ = GHZ*GHZ';
% W = W*W';
% rho = 0.6*GHZ + 0.4*W;
% ge_mixed_ppt_kext(rho, [2,2,2], 2)   % example with 2-symmetric extension
% ans =
%    0.30252001211618484
    

    if nargin < 4 || isempty(sdpaccuracy)
        sdpaccuracy = 1e-8;
    end
    if nargin < 5 || isempty(itera)
        itera = 200;
    end

    n=length(dim);
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
    if floor(k) ~= k
        error('k must be an integer');
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
    if n==2              %k extendibility for bipartite states
        d0 = dim(1);
        d1 = dim(2);
        pro = d0 * d1;
        
        
        
        if k<=1
            error('k must be at least 2')
        end
        
        n_dim = d0 * d1^k;
        
        
        cvx_begin sdp
            cvx_solver sdpt3
            cvx_solver_settings( 'gaptol', sdpaccuracy, 'inftol', sdpaccuracy,  'steptol', sdpaccuracy,  'maxit', itera)
    
            variable siigma(n_dim, n_dim) hermitian semidefinite
            variable X(pro, pro) complex;
            
            
            
            sigpt=PartialTrace(siigma,3:k+1,[d0, repmat(d1, 1, k)])
            
            for i=2:k
                
                wa=1:k;
                wa(i)=[];
                PartialTrace(siigma,wa+1,[d0, repmat(d1, 1, k)])==sigpt;
            end
            trace(sigpt) == 1;
            sigpt >= 0;
            for i=1:k+1
                PartialTranspose(siigma, i, [d0, repmat(d1, 1, k)])>=0; 
                
            end
            
            B = [rho, X; X', sigpt];
            B >= 0;
            
          
            A = [zeros(pro), eye(pro)/2; eye(pro)/2, zeros(pro)];
            maximize trace(A * B)
        cvx_end
        
       
    if strcmp(cvx_status, 'Inaccurate/Solved')
        warning('Solution did not match with desired precision. Solution may be inaccurate.');
        
        result = [1 - (cvx_optval+sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)))^2,0];
    else
        result = [1 - (cvx_optval+sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)))^2,1];
    end
    end
    if n>2          %k extendibility for multipartite states
        reszta = mod(k, n - 1);
        kk = floor(k / (n - 1));
        
        ran = 0:(length(dim)-1);
        ranco = ran;
        ranco = fliplr(ranco);
        n_dim = [dim, repmat(dim(2:end), 1, kk), dim(n - reszta + 1:end)];
        n_num = [ran, repmat(ran(2:end), 1, kk), ran(n - reszta + 1:end)];
        pro=prod(dim);
        n_pro=prod(n_dim);
        n_number=length(n_dim);
        add_di = ones(1, n-1);
        for i = 0:(k-1)
            add_di(mod(i, n-1) + 1) = add_di(mod(i, n-1) + 1) + 1;
        end
        
        cvx_begin sdp
            cvx_solver sdpt3
            cvx_solver_settings( 'gaptol', sdpaccuracy, 'inftol', sdpaccuracy,  'steptol', sdpaccuracy,  'maxit', itera)
    
            variable siigma(n_pro, n_pro) hermitian semidefinite
            variable X(pro, pro) complex;
            for sy = 1:(n-1)   
                
                for napot = 1:(add_di(sy) - 1)
                    ile = 0;
                    ilesys = 0;
                    zostaw = [];
                    zostaw2 = [];
                    czy = 1;
                    wywal = [];
                    
                    % while len(zostaw+zostaw2) < n
                    while length([zostaw, zostaw2]) < n
                        ile = ile+1;
                        
                        
                        if n_num(ile) ~= ranco(sy) && czy == 1
                            zostaw = [zostaw, ile];
                        
                        elseif n_num(ile) ~= ranco(sy) && czy == 0
                            wywal = [wywal, ile];
                        
                        else
                            ilesys = ilesys + 1;
                            
                            if ilesys <= napot
                                czy = 0;
                                wywal = [wywal, ile];
                            end
                            
                            if ilesys == napot + 1
                                zostaw2 = [zostaw2, ile];
                                
                                if sy > 1
                                    for isy = 1:sy-1
                                        zostaw2 = [zostaw2, ile + isy];
                                    end
                                end
                            end
                        end
                        
                    end
                    n_wywal=1:n_number;
                    
                    n_wywal(unique([zostaw, zostaw2])) = [];
                    indeksy=1:length(n_dim);
                    
                    indeksy(n+1:end);
                    PartialTrace(siigma,n_wywal,n_dim)==PartialTrace(siigma,indeksy(n+1:end),n_dim);

                end
            end
        
                        
            sigpt=PartialTrace(siigma,n+1:n_number,n_dim);
            
            
            trace(sigpt) == 1;
            sigpt >= 0;
            for i=1:n_number
                PartialTranspose(siigma, i, n_dim)>=0; 
                
                
            end
            
            B = [rho, X; X', sigpt];
            B >= 0;
            
          
            A = [zeros(pro), eye(pro)/2; eye(pro)/2, zeros(pro)];
            maximize trace(A * B)
        cvx_end
        
    if strcmp(cvx_status, 'Inaccurate/Solved')
        warning('Solution did not match with desired precision. Solution may be inaccurate.');
        
        result = [1 - (cvx_optval+sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)))^2,0];
    else
        result = [1 - (cvx_optval+sdpaccuracy*(1+2*(cvx_optval+sdpaccuracy)))^2,1];
    end
        
    end

end