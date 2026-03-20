function Op = TensorWithout(phiCell, idx)
% TensorWithIdentityAt  Constructs tensor product from list of kets,
% replacing the idx-th ket with an identity operator.
%
% INPUT:
%   phiCell : 1xN cell array, each entry is a ket vector (column vector)
%   idx     : index (1-based) of ket to replace with identity
%
% OUTPUT:
%   Op      : operator (matrix), tensor product with identity at idx

    N = length(phiCell);
    factors = cell(1,N);
    for j = 1:N
        if j == idx
            d = length(phiCell{j});        % dimension of subsystem j
            factors{j} = eye(d);           % identity instead of ket
        else
            v = phiCell{j} / norm(phiCell{j});  % normalize ket
            factors{j} = v ;           % projector |v><v|
        end
    end
    Op = Tensor(factors);  % QETLAB Tensor
end