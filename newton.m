% Newton's method to find x s.t. F(x) = 0
function [xk,iter, err] = newton (F,x0, x, eps)

xk = x0;
Fk = subs(F,x,xk);
jac = jacobian(F, x);
Jk = subs(jac, x, xk);
ek = double(norm(Fk));
iter = 0;
err = ek;
while(ek > eps)
    
    xk = xk-(Jk\Fk);
    Fk = double(subs(F, x, xk));
    Jk = subs(jac, x, xk);
    ek = double(norm(Fk));
    err = [err, ek];
    iter = iter + 1;
    if iter > 200
        break
    end
    
end


end