% Broyden's method to find x s.t. F(x) = 0
function [xk,iter, err] = broyden(F,x0, A0, x, eps)

xk = x0;
Ak = A0;
Fk = subs(F, x, xk);
ek = double(norm(Fk));
iter = 0;
err = ek;
while (ek > eps)
    
    sk = -Ak\Fk; % solve the linear system Ak*sk = -Fk
    xk = xk +sk; % update xk value
    yk = subs(F, x, xk) - Fk;
    Fk = subs(F, x, xk); % update Fk value
    Ak = Ak + 1/(transpose(sk)*sk)*(yk-Ak*sk)*transpose(sk); % Broyden's update
    ek = double(norm(Fk));
    err = [err, ek];
    iter = iter + 1;
    if iter > 200
        break
    end
end

end