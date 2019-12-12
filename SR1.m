% Symetric rank 1 quasi-Newton method
function [xk,iter,err] = SR1 (f,gradf,H0, x, x0, eps)
    
xk = x0;
Hk = H0;
Bk = H0^(-1) ;
iter = 0;
gk = double(subs(gradf,x,xk));
ek = norm(gk);
err = [];
while(ek>eps)
    sk     = -Bk*gk;
    xk     = sk' + xk;
    gk_new = double(subs(gradf,x,xk)); 
    ek     = norm(gk_new);
    err    = [err,ek];
    yk     = gk_new - gk;
    Bk     = Bk + (sk-Bk*yk)*(sk-Bk*yk)'*(1/((sk-Bk*yk)'*yk));
    iter   = iter + 1;
    gk     = gk_new;
    if iter > 200
        break
    end
end

end
