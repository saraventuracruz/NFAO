% Broyden's method to find x s.t. F(x) = 0
function [xk,iter, err] = broyden(F,x0, A0, x, eps)

xk = x0; % initialize xk
Ak = A0; % initialize Ak
Fk = subs(F, x, xk); % compute F0
ek = double(norm(Fk)); % compute current error
iter = 0; % initalize iterations' counter
err = ek; % save the current error for output
while (ek > eps)
    
    sk = -Ak\Fk; % solve the linear system Ak*sk = -Fk
    xk = xk +sk; % update xk value
    yk = subs(F, x, xk) - Fk;
    Fk = subs(F, x, xk); % update Fk value
    Ak = Ak + 1/(transpose(sk)*sk)*(yk-Ak*sk)*transpose(sk); % Broyden's update
    ek = double(norm(Fk)); % compute current error
    err = [err, ek]; % save the current error for output
    iter = iter + 1; % update the iterations' counter
    if iter > 200
        break
    end
end

end