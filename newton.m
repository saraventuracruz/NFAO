% Newton's method to find x s.t. F(x) = 0
function [xk,iter, err] = newton (F,x0, x, eps)

xk = x0; % initialize xk
jac = jacobian(F, x); % define the jacobian of F
Fk = double(subs(F, x, xk)); % compute the first Fk
Jk = subs(jac, x, xk); % evaluate the jacobian on x0
ek = double(norm(Fk)); % update the error
err = ek; % % save the first error for output
iter = 0; % initialize iterations' counter
while(ek > eps)
    
    
    
    xk = xk-(Jk\Fk); % update xk
    Fk = double(subs(F, x, xk)); % compute current Fk
    Jk = subs(jac, x, xk); % compute current jacobian
    ek = double(norm(Fk)); % update the error
    err = [err, ek]; % save the current error for output
    iter = iter + 1; % update the iterations' counter
    if iter > 200
        break
    end
    
end


end