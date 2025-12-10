ghz_state = [1; 0; 0; 0; 0; 0; 0; 1]; 
ghz_state = ghz_state / norm(ghz_state);


q5=RandomStateVector(32);

%The second position at the output, the error of relaxing separability
%condition to PPT is 0. It is due to the fact that we optimized over
%two-qubit system and in this case set of separable states is equal
%to set of PPT states. For higher dimesnions there is an
%estimation error
%caused by difference between separable and PPT set of states.

res=ge_pure(ghz_state,[2,2,2]);
res

res5q=ge_pure(q5,[2,2,2,2,2]);
res5q