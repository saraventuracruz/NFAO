function [XK, xk,lambdak,iter, cond_H] = augmented_lagrangian_method_Newton(x, x0, objFunction, constraints, mu0,lambda0, eps)

xk = transpose(x0);
muk = mu0;
lambdak = lambda0;
cond_H = [];
iter = 0;
XK = [];
while(true)
    iter = iter +1;
    AL = augmented_lagrangian(lambdak, objFunction, constraints, muk);
    grad_al = gradient(AL, x); % gradient
    hess_al = hessian(AL, x); % hessian
   
    H0_al = double(subs(hess_al,transpose(x),xk));
    cond_H = [cond_H;real(double(cond(H0_al)))];
    [~,xk,~,~] = pure_newton(AL,grad_al,hess_al,transpose(x), xk, eps);
    XK = [XK;xk];
    t = 0;
    for i = 1:size(constraints, 1)
        t = t+abs(double(subs(constraints(i), x, transpose(xk))));
    end
    if(t <= eps)
        break;
    end
    lambdak = lambdak-double(subs(constraints, x, transpose(xk)))/muk;
    muk = muk*0.1;

end
end