function [xk,iter, err] = newton (F,x0, var_x, var_lambda, eps)

xk = x0;
Fk = subs(F,[transpose(var_x) transpose(var_lambda)],xk)
jac = jacobian(F, [transpose(var_x) transpose(var_lambda)]);
Jk = subs(jac, [transpose(var_x) transpose(var_lambda)], xk)
ek = norm(Fk);
iter = 0;
err = ek;
while(ek > eps)
    
    xk = xk-(Jk\Fk)';
    Fk = double(subs(F, [transpose(var_x) transpose(var_lambda)], xk));
    Jk = subs(jac, [transpose(var_x) transpose(var_lambda)], xk);
    ek = norm(Fk);
    err = [err, ek];
    iter = iter + 1;
    
end


end