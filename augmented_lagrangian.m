function l = augmented_lagrangian(lambda, objFunction, constraints, mu)
    
    l = objFunction;
    for i = 1:size(constraints)
        l = l - lambda(i)*constraints(i) + 1/(2*mu)*constraints(i)^2;
    end

end