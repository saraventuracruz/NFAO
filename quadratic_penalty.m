function q = quadratic_penalty(objFunction, constraints, mu)

    q = objFunction;
    for i = 1:size(constraints)
        q = q+1/(2*mu)*constraints(i)^2;
    end

end