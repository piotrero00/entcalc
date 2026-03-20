function output = upperbip(rho,dims,dec,insep,r,iteramax,dif)
% upperbip  Computes the upper bound of the geometric entanglement 
% of a bipartite state rho.
%
% Parameters
% ----------
% rho : matrix
%     A 2D array representing a mixed state.
%
% dims : vector of integers
%     Dimensions of each subsystem. Each entry should be a positive
%     integer. Number of elements in dims must be 2.
%
% dec : integer (0 or 1)
%     If 1, the separable decomposition is also returned.
%     If 0, only the bound is returned.
%
% insep : cell
%     Cell array that specifies the initial separable decomposition.
%     - If the first element is 0, the cell contains only this element
%       and the initial separable decomposition will be chosen randomly.
%     - If the first element is 1, the cell must contain three elements:
%         {1, ps, kets}
%       where:
%         ps   : vector of probabilities corresponding to the kets
%         kets : cell array of pure states (column vectors) used in the
%                initial separable decomposition
%
% r : integer, optional
%     Number of pure states in the separable decomposition of rho.
%     Default is prod(dim)^2, where dim are the subsystem dimensions of rho.
%     A recommended choice is (rank(rho))^2.
%
% iteramax : integer, optional
%     Maximum number of iterations of the algorithm. Default is 1000.
%
% dif : float, optional
%     Convergence threshold for the algorithm. Default is 1e-8.
%
% Returns
% -------
% output : float or cell
%     If dec = 0 (default), the function returns a single float value –
%     the upper bound on the geometric entanglement of rho.
%
%     If dec = 1, the function returns a cell array with three elements:
%         output{1} – upper bound (float)
%         output{2} – vector of probabilities corresponding to the
%                     separable decomposition
%         output{3} – cell array of pure states (column vectors)
%
%
% Notes
% -----
% - The function computes the upper bound using an iterative algorithm.
%   Each iteration refines the separable decomposition of the maximization variable.
%
% - The input state must be bipartite.
%
% - Due to square roots and machine precision (~1e-16), the numerical accuracy
%   of the computation is around 1e-8. This tolerance is not explicitly
%   included in the result, meaning that the returned upper bound can be up
%   to 1e-8 lower than the true geometric entanglement.
%
% Examples
% --------
% GHZ = [1; 0; 0; 0; 0; 0; 0; 1] / sqrt(2);
% W = [0; 1; 1; 0; 1; 0; 0; 0] / sqrt(3);
% GHZ = GHZ*GHZ';
% W = W*W';
% rho = 0.6*GHZ + 0.4*W;
% upperbip(rho,[2,2,2],0,{0})
% ans =
%    0.2892



    n=prod(dims);
    if nargin < 5 || isempty(r)
        r = n^(2);
    end
    if nargin < 6 || isempty(iteramax)
        iteramax = 1000;
    end
    if nargin < 7 || isempty(dif)
        dif = 10^(-8);
    end
    pqp=eig(rho);
    licznik = sum(pqp > 1e-15);
    
    if abs(trace(rho) - 1) > 1e-14
        error('State should be normalized');
    end
    
    
    if length(dims) ~= 2
        error('System should be a bipartite quantum state');
    end
    
    
    if min(dims) <= 1
        error('Dimension of subsystem must be an integer greater than 1');
    end
    
    
    if any(floor(dims) ~= dims)
        error('Dimension of subsystem must be an integer greater than 1');
    end
    
    
    if norm(rho - rho', 'fro') > 1e-12
        error('State must be represented by a positive-semidefinite matrix.');
    end
    
    
    evals = eig(rho);
    if any(evals < -1e-10)
        error('State must be represented by a positive-semidefinite matrix.');
    end
    
    
    if size(rho,1) ~= dims(1)*dims(2)
        error('Dimension of the system is not equal to the product of dimensions of subsystem');
    end
    if r> (dims(1)*dims(2))^(2)
        error('r should be maximally dim(rho)^(2)');
    end
    if r<licznik
        error("r is smaller than rank(rho). r should be at least equal to rank(rho)")
    end
    if r<licznik^(2)
        warning("r is smaller than rank(rho)^(2). The output might be less precise")
    end
    if ~(dec==0 || dec==1)
        error('dec must be 0 or 1');
    end
    if insep{1}==1
        if length(insep{2})~=r
            error('r should have the same number of elements as sep{2}')
        end
        if length(insep{3})~=r
            error('r should have the same number of elements as sep{3}')
        end
    end
    if floor(r)~=r
        error('r should be an integer')
    end
    if floor(iteramax) ~= iteramax
        error('iteramax must be an integer greater than 0');
    end
    if iteramax <= 0
        error('iteramax must be an integer greater than 0');
    end
    if ~isa(dif,'double')
        error('dif must be a float greater than 0');
    end
    if dif >=1 
        error('dif must be a float smaller than 1');
    end
    if dif <= 0
        error('dif must be a float greater than 0');
    end

    %Main part of the function
    p = cell(r,1);
    psi = cell(r,1);
    
    
    
    e = eye(r);
    
    %Creating a random decomposition of rho
    
    U = RandomUnitary(r);
    
    [V,D] = eigs(rho,r);
    
    
    %{p,psi} is a pure state decomposition of the input state rho
    for i = 1:r
        j = 1;
        psi{i} = U(i,j)*sqrt(D(j,j))*V(:,j);
    
        for j=2:min(r,prod(dims))
            psi{i} = psi{i} + U(i,j)*sqrt(D(j,j))*V(:,j);
        end
        p{i} = psi{i}'*psi{i};
        psi{i} = psi{i}/sqrt(p{i});
    
    end
    

    %{q,phi} will be a separable decomposition of the closest separable state
    if insep{1}==0
        q = rand(r,1);
        q = q/sum(q);
        q = num2cell(q);
        phiCell = cell(r,1);
        rvCell = cell(length(dims),1);
        
        
        for i = 1:r
            for j = 1:length(dims)
                rvCell{j} = RandomStateVector(dims(j));
            end
        
            phiCell{i} = rvCell;
        
        end
    else
        q=insep{2};
        phiCell=insep{3};
    end
    F=42;
    Fp=0;
    %Iterations, which improves separable decomposition
    for l = 1:iteramax
        if abs(F-Fp)<dif
            break
        end
        A = zeros(r);
    
        for i = 1:r
            for j = 1:r
                A = A + sqrt(p{i}*q{j}) * Tensor(phiCell{j})'*psi{i} * e(:,i)*e(:,j)';
            end
        end
    
        [V,~,W]=svd(A);
        U = W*V';
    
        alpha = cell(r,1);
        for i = 1:r
            alpha{i} = zeros(prod(dims),1);
            for j = 1:r
                alpha{i} = alpha{i} + U(i,j)*sqrt(p{j})*psi{j};
            end
        end
    
        pNew = cell(r,1);
        psiNew = cell(r,1);
    
        for i = 1:r
            pNew{i} = alpha{i}'*alpha{i};
            psiNew{i} = alpha{i}/sqrt(pNew{i});
        end
    
        phiNewCell = cell(r,1);
        qNew = cell(r,1);
    
        for i = 1:r
            phiNewCell{i} = bipup(psiNew{i},phiCell{i},dims);  %Improving separable decomposition
        end
    
        
        summ = 0;
        for k = 1:r
            summ = summ + pNew{k}*abs(psiNew{k}'*Tensor(phiNewCell{k}))^2;
        end
    
        for i = 1:r
            qNew{i} = pNew{i} * abs(psiNew{i}'*Tensor(phiNewCell{i}))^2/summ;
        end
    
    
    
        psi = psiNew; phiCell = phiNewCell; p = pNew; q = qNew;
        
        si=0;
        for iu=1:r
            si=si+qNew{iu}*Tensor(phiNewCell{iu})*Tensor(phiNewCell{iu})';
        end
        rin=sqrtm(si)*rho*sqrtm(si);
    
        rsq=sqrtm(rin);
        Fp=F;
        F=1-real(trace(rsq)^(2));      %Computing fidelity

    end
        if dec==1
            output=cell(3,1);
            
            output{1}=F;
            output{2}=q;
            output{3}=phiCell;
        else
            output=F;
    
        end
    
    end
    

