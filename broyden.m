function [Xk, xk,iter,err] = broyden (F,A0,x0, eps)

Ak = A0;
xk = x0;
while

sk = -inv(A0)*F(xk);


yk = F(xk+sk)-F(xk);
xk = xk + sk;

Ak = Ak + 1/(sk'*sk)*(yk-Ak*sk)*sk';

end

end