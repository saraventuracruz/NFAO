function l = lagrangian (lambda, mu, objFunction, eqConst, ineqConst)
    
    nEqConst = size(eqConst);
    nIneqConst = size(ineqConst);
    
    l = objFunction;
    for i = 1:nEqConst
        l = l - lambda(i)*eqConst(i);
    end
    for i = 1:nIneqConst
        l = l-mu(i)*ineqConst(i);
    end
end