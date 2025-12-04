function mpCell = bipup(psi,phiCell,dims)

    mpCell = phiCell;



    theta0 = PartialTrace(psi*psi', 2, dims);  
    
    
    theta1 = PartialTrace(psi*psi', 1, dims);  
    
   
    [V0,D0] = eig(theta0); 
    [V1,D1] = eig(theta1); 
    
    
    [~,idx0] = max(diag(D0));
    [~,idx1] = max(diag(D1));
    

    cc = V0(:,idx0);
    dd = V1(:,idx1);
    mpCell = {cc, dd};


end