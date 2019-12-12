% NFAO 
% Computational Assignment - Version 5
% Kilian & Sara

% 1) Unconstrained optimization 
%% Question a : Prove that x_star is stationary
N   = 4; % dimension of the problem
x = sym('x', [1,N]) % definition of the variable x (row vector)
f = fun(x) % f is the function defined on fun.m
grad_f = gradient(f,x) % compute the gradient of f using the matlab function gradient
x_star = [1, 1, 1, 1]; % row vector
f_star = subs(f,x,x_star) % f applied to x_star (= 0)
grad_star = subs(grad_f,x,x_star) % gradf applied to x_star (= [0,0,0,0])
% grad(f) vanishes in x_star => x_star stationary

%% Question b : Prove necessary first and second order optimality conditions
hess_f = hessian(f,x) % compute the hessian matrix of f using the matlab function
hess_star = double(subs(hess_f,x,x_star)) % hessian apllied to x_star
[~,D] = eig(hess_star); 
eigenvalues = sort(diag(D)) % eigenvalues are non negative
% first order : x_hess_star = double(subs(hessf,x,x_star))
[~,D] = eig(hess_star);
eigenvalues = sort(diag(D)) % star stationary
% second order : hess(f) in x_star is symetric defined positive

%% Question c : 
% f(x) >= 0 for all x and f(x) = 0 iff x = x_star
kkt_system = grad_f == zeros(N,1)
x_star = vpasolve(kkt_system, x); % solve numerically
% get values 
fieldNames = fieldnames(x_star); %x1 x2 ...
sol = []; % to store each solution
for i= 1:size(fieldNames)
    sol(i,:) = x_star.(fieldNames{i})'; % access the solution and save it transposed (in order to get a matrix such that the number of lines is the problem dimension and the number of columns is the number of different solutions)
end
% get real solutions only
stationary_points = [];
for i = 1:size(sol, 2)
    if isreal(sol(:,i))
        stationary_points = [stationary_points, sol(:, i)];
    end
end
stationary_points % print the stationary points to screen

%% Question d :

x0 = rand(1, N); % random row vector of dimension N
eps = 1e-6; % epsilon
tic()
[X_newton, x_newton,iter_newton,err_newton] = pure_newton (f,grad_f,hess_f, x, x0, eps);
x_newton, iter_newton
grad_newton = double(subs(grad_f,x,x_newton))
elapsed_time_d = toc()

%% Question f :

tic()
x0 = rand(1, N);
gamma = 0;
H0 = double(subs(hess_f,x,x0)) + gamma*eye(N);
[x_SR1,iter_SR1,err_SR1] = SR1 (f,grad_f,H0, x, x0, eps);
x_SR1, iter_SR1
grad_SR1 = double(subs(grad_f,x,x_SR1))
elapsed_time_f = toc()

%% Question g : 
figure(1)
xx = 1:iter_newton
yy = err_newton
plot(xx,yy,'x')
title('||graf(x_k)|| newton')

figure(2)
xx = 1:iter_SR1
yy = err_SR1
plot(xx,yy,'x')
title('||graf(x_k)|| SR1')


%%
% 2) Constrained Optimization
%
%% Question a :
%% Problem definition

N = 3; %problem dimension
x = sym('x',  [N,1]); % definition of variable x (column vector of dimension N)
objFun = fun2(x) % definition of the objective function (using the function defined on fun2.m)
eqConst = [cons1(x); cons2(x)] % definition of equality constraints (using the functions defined on cons1.m and cons2.m)
ineqConst = []; % definition of inequality constraints

%% Definition of the KKT system

nEqConst = size(eqConst,1); % number of equality constraints
nIneqConst = size(ineqConst,1); % number of inequality constraints
lambda = sym('lambda', [nEqConst,1]); % definition of variable lambda (set of lagrangian multipliers for the equality constraints)
mu = sym('mu', [nIneqConst,1]); % definition of variable mu (set of lagrangian multipliers for the inequality constraints)
lag = lagrangian(lambda, mu, objFun,eqConst,ineqConst) % definition of the lagrangian function (using the function defined at lagrangian.m)
grad_lag = gradient(lag, x); % definition of the gradient of the lagrangian function

kkt_system = [simplify(grad_lag); eqConst; ineqConst]; % definition of the system of KKT system
% for i = 1:nEqConst
%     kkt_system = [kkt_system; lambda(i)*eqConst(i)];    
% end

%nEquations = N + 2*nEqConst; % Total number of equations on the KKT system
nEquations = N + nEqConst;
kkt_conditions = kkt_system == zeros(nEquations, 1) % print the equations to screen

% inequality constraints conditions:
for i = 1:nIneqConst   
    kkt_conditions = [kkt_conditions; ineqConst(i) >= 0];
    kkt_conditions = [kkt_conditions; mu(i)*ineqConst(i) == 0];
    kkt_conditions = [kkt_conditions; mu(i)>=0]; 
end

%% Solve the KKT system
KKT_sol = solve(kkt_conditions, [transpose(x) transpose(lambda) transpose(mu)])
% get values 
fieldNames = fieldnames(KKT_sol); %x1 x2 ...
x_star = []; % to store each solution
for i= 1:N
    x_star(i,:) = KKT_sol.(fieldNames{i})'; % access the solution and save it transposed (in order to get a matrix such that the number of lines is the problem dimension and the number of columns is the number of different solutions)
end
lambda_star = [];
for i = N+1:nEqConst+N
    lambda_star(i-N, :) = KKT_sol.(fieldNames{i})';
end
x_star,lambda_star


%% Question d:
y0 = rand(1,nEquations);
[xk_newton, iter_newton, err_newton] = newton(kkt_system, y0, x, lambda, 1e-6)

%% Question e:
mu = 0.5;
q = quadratic_penalty(objFun, eqConst, mu)
grad_q = gradient(q, x)
hess_q = hessian(q, x)

q0 = rand(1, N);
pure_newton(q, grad_q, hess_q, x, q0, eps)

[Q_newton, q_newton,q_iter_newton,q_err_newton] = pure_newton (q,grad_q,hess_q, x, q0, eps);
q_newton, q_iter_newton

gamma = 0;
H0_q = double(subs(hess_q,x,q0)) + gamma*eye(N)
[q_x_SR1,q_iter_SR1,q_err_SR1] = SR1 (q,grad_q,H0_q, x, q0, eps);
q_x_SR1, q_iter_SR1

