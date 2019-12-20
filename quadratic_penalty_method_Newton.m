function [XK,xk, iter, cond_H] = quadratic_penalty_method_Newton(x, x0, objFunction, constraints, mu0, eps)

xk = transpose(x0);
muk = mu0;
cond_H = [];
iter = 0;
XK = [];
while(true)
    iter = iter +1;
    Q = quadratic_penalty(objFunction, constraints, muk);
    grad_q = gradient(Q, x); % gradient
    hess_q = hessian(Q, x); % hessian
        
    H0_q = double(subs(hess_q,transpose(x),xk));
    cond_H = [cond_H;real(double(cond(H0_q)))];
    [~,xk,~,~] = pure_newton(Q,grad_q,hess_q,transpose(x), xk, eps);
    XK = [XK;xk];
    t = 0;
    for i = 1:size(constraints, 1)
        t = t+1/(2*muk)*double(subs(constraints(i), x, transpose(xk)))^2;
    end
    if(t <= eps)
        break;
    end
    muk = muk*0.1;
end
end