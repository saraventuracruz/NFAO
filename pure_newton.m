%solves gradf(x) = 0
function [Xk, xk,iter,err] = pure_newton (f,gradf,hessf, x, x0, eps)
    
xk = x0;
% fk = double(subs(f,x,xk));
gk = double(subs(gradf,x,xk));
ek = norm(gk);
Hk = double(subs(hessf,x,xk));
iter = 0;
err = [];
Xk = x0;
while(ek>eps)
    sk = -Hk\gk;
    xk = sk' + xk;
    Xk = [Xk;xk];
    gk_new = double(subs(gradf,x,xk));
    ek = norm(gk_new);
    err = [err,ek];
    gk = gk_new;
    Hk = double(subs(hessf,x,xk));
    iter = iter + 1;
end


end
