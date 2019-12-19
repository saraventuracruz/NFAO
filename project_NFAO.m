% NFAO 
% Computational Assignment - Version 5
% Kilian & Sara

%% define 
format long
%% 1) Unconstrained optimization 
%% Question a : Prove that x_star is stationary
N   = 4; % dimension of the problem
x = sym('x', [1,N]) % definition of the variable x (row vector)
f = fun1(x) % f is the function defined on fun1.m
grad_f = gradient(f,x) % compute the gradient of f using the matlab function gradient
x_star = [1, 1, 1, 1]; % row vector
f_star = subs(f,x,x_star) % f applied to x_star (= 0)
grad_star = subs(grad_f,x,x_star) % gradf applied to x_star (= [0,0,0,0])
% grad(f) vanishes in x_star => x_star stationary

%% Question b : Prove necessary first and second order optimality conditions
hess_f = hessian(f,x) % compute the hessian matrix of f using the matlab function
hess_star = double(subs(hess_f,x,x_star)) % hessian apllied to x_star
eigenvalues = eig(hess_star); % eigenvalues are non negative
% second order : hess(f) in x_star is symetric defined positive

%% Question c : 
% f(x) >= 0 for all x and f(x) = 0 iff x = x_star
eq = grad_f == zeros(N,1)
x_star = vpasolve(eq, x); % solve numerically
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

% Compute the eigenvalues of the hessian apllied on the stationary points
eig_values = [];
for i = 1:size(stationary_points, 2)
    eig_values = [eig_values, eig(double(subs(hess_f,x,stationary_points(:,i)')))];
end

eig_values
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
[X_SR1,x_SR1,iter_SR1,err_SR1] = SR1 (grad_f,H0, x, x0, eps);
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
%% Problem definition

N = 3; % problem's dimension
x = sym('x',  [N,1]); % definition of variable x (column vector of dimension N)
objFun = fun2(x) % definition of the objective function (using the function defined on fun2.m)
eqConst = [cons1(x); cons2(x)] % definition of equality constraints (using the functions defined on cons1.m and cons2.m)
ineqConst = []; % definition of inequality constraints

%% Test definiteness of A
A = [ 6, 2, 1; 2, 5, 2; 1, 2, 4]; % define matrix A
eig_A = eig(A); % compute eigenvalues of A
eig_A > 0 % test positive sign of each eigenvalue

%% Question a:
%% Definition of the KKT system

nEqConst = size(eqConst,1); % number of equality constraints
nIneqConst = size(ineqConst,1); % number of inequality constraints
lambda = sym('lambda', [nEqConst,1]); % definition of variable lambda (set of lagrangian multipliers for the equality constraints)
mu = sym('mu', [nIneqConst,1]); % definition of variable mu (set of lagrangian multipliers for the inequality constraints)
lag = lagrangian(lambda, mu, objFun,eqConst,ineqConst) % definition of the lagrangian function (using the function defined at lagrangian.m)
grad_lag = gradient(lag, x); % definition of the gradient of the lagrangian function

% definition of the system of KKT system
kkt_system = [simplify(grad_lag); eqConst;ineqConst]; 
for i = 1:nEqConst
    kkt_system = [kkt_system; lambda(i)*eqConst(i)];    
end

nEquations = N + 2*nEqConst; % Total number of equality conditions on the KKT system
kkt_conditions = kkt_system == zeros(nEquations, 1); % definition of the kkt conditions (= 0)

% add the inequality constraints conditions, when there are inequality constraints:
for i = 1:nIneqConst   
    kkt_conditions = [kkt_conditions; ineqConst(i) >= 0];
    kkt_conditions = [kkt_conditions; mu(i)*ineqConst(i) == 0];
    kkt_conditions = [kkt_conditions; mu(i)>=0]; 
end

kkt_conditions % print out the KKT conditions

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
% redefine kkt_system for the current problem (drop out the conditions of
% lambda * eqConst = 0
kkt_system = kkt_system(1:(end-nEqConst),1)
%% Solve using Newton's method
z0 = randi(10, size(x,1) + size(lambda,1), 1) %random initial guess (entries are integers between 0 and 10)
eps = 1e-6; %epsilon
[zk_newton, iter_newton, err_newton] = newton(kkt_system, z0, [x; lambda], eps)
% xk_newton is the concatenation of x_star and lambda_star found by
% newton's method
%% Solve using Broyden's method
z0Broyden = [];
iterBroyden = [];
zStarBroyden = [];
for i = 1:5
    z0 = randi(10,size(x,1) + size(lambda,1), 1); %random initial guess (entries are integers between 1 and 10)
    z0Broyden = [z0Broyden; transpose(z0)];
    eps = 1e-8; %epsilon
    A0 = rand(size(x,1) + size(lambda,1)); % generate a random squared matrix with the same dimension of the kkt system
    A0 = A0*transpose(A0); % make it positive semidefinite
    [zk_broyden, iter_broyden, err_broyden] = broyden(kkt_system, z0, A0, [x; lambda], eps);
    iterBroyden = [iterBroyden; iter_broyden];
    zStarBroyden = [zStarBroyden; transpose(zk_broyden)];
end

%% Print out table with results
T_Broyden = table(z0Broyden, iterBroyden, zStarBroyden)


%% Question e:
mu = 0.5;
Q = quadratic_penalty(objFun, eqConst, mu) % definition of the quadratic penalty method (quadratic_penalty.m)
grad_q = gradient(Q, x) % gradient
hess_q = hessian(Q, x) % hessian


%% Solve using Pure Newton's method
q0 = rand(1, N)
eps = 1e-6;
[Q_newton, q_newton,q_iter_newton,q_err_newton] = pure_newton (Q,grad_q,hess_q, transpose(x), q0, eps);
q_newton, q_iter_newton

%% Solve using SR1
q0 = rand(1, N)
eps = 1e-6;
gamma = 5;
H0_q = double(subs(hess_q,transpose(x),q0)) + gamma*eye(N)
[~,q_x_SR1,q_iter_SR1,q_err_SR1] = SR1 ( grad_q, H0_q, transpose(x), q0, eps);
q_x_SR1, q_iter_SR1

%% test different mu values
q0 = rand(1, N);
%% Solve using Pure Newton's method
results_Q_newton = [];
for i = -4:10
    mu = 10^i;
    Q = quadratic_penalty(objFun, eqConst, mu); % define the quadratic penalty function for the current mu
    grad_q = gradient(Q, x); % gradient
    hess_q = hessian(Q, x); % hessian
    % Solve using pure newton's method
    [~, q_newton, q_iter_newton, q_err_newton] = pure_newton(Q, grad_q, hess_q, transpose(x), q0, eps);
    results_Q_newton = [results_Q_newton; [q_newton,q_iter_newton, q_err_newton(end)]];
end
%% Print out table with results
mu = transpose(10.^(-4:10));
xStarNewton = results_Q_newton(:, 1:N);
iterNewton = results_Q_newton(:, N+1);
errorNewton = results_Q_newton(:, N+2);
T_Q_Newton = table(mu, xStarNewton, iterNewton, errorNewton)
%% Solve using SR1
results_Q_SR1 = [];
for i = -4:2:10
    mu = 10^i;
    Q = quadratic_penalty(objFun, eqConst, mu); % define the quadratic penalty function for the current mu
    grad_q = gradient(Q, x); % gradient
    hess_q = hessian(Q, x); % hessian
    for gamma = 0:5:10
        H0_q = double(subs(hess_q,transpose(x),q0)) + gamma*eye(N);
        [Q_x_SR1,q_x_SR1,q_iter_SR1,q_err_SR1] = SR1 (grad_q, H0_q, transpose(x), q0, eps);
        results_Q_SR1 = [results_Q_SR1; [q_x_SR1, q_iter_SR1,q_err_SR1(end)]];
    end
end
%% Print out table
mu = [];
gamma = [];
for i = -4:2:10
    for g = 0:5:10
        mu = [mu; 10^i];
        gamma = [gamma ; g];
    end
end
xStarSR1 = results_Q_SR1(:, 1:N);
iterSR1 = results_Q_SR1(:, N+1);
errorSR1 = results_Q_SR1(:, N+2);
T_Q_SR1 = table(mu, gamma, xStarSR1, iterSR1, errorSR1)

%% Question f:
mu = 0.2;
ALag = augmented_lagrangian(lambda, objFun, eqConst, mu) %definition of the augmented lagrangian method (augmented_lagrangian.m)
grad_aLag = gradient(ALag, [transpose(x) transpose(lambda)]) % gradient
hess_aLag = hessian(ALag, [transpose(x) transpose(lambda)]) % hessian

%% Solve using Pure Newton's method
q0 = rand(1, N+nEqConst);
eps = 1e-6;
[ALag_newton, aLag_newton, aLag_iter_newton, aLag_err_newton] = pure_newton (ALag,grad_aLag,hess_aLag, [transpose(x) transpose(lambda)], q0, eps);
aLag_newton, aLag_iter_newton

%% Solve using SR1
q0 = rand(1, N+nEqConst);
eps = 1e-6;
gamma = 0;
H0_aLag = double(subs(hess_aLag,[transpose(x) transpose(lambda)],q0)) + gamma*eye(N+nEqConst)
[aLag_x_SR1,aLag_iter_SR1,aLag_err_SR1] = SR1 (grad_aLag,H0_aLag, [transpose(x) transpose(lambda)], q0, eps);
q_x_SR1, q_iter_SR1

%% test different mu values
l0 = rand(1, N+nEqConst);

%% Solve using Pure Newton's method
results_L_newton = [];
for i = -4:10
    mu = 10^i;
    ALag = augmented_lagrangian(lambda, objFun, eqConst, mu); % define the augmented lagrangian for the current mu
    grad_aLag = gradient(ALag, [transpose(x) transpose(lambda)]); % gradient
    hess_aLag = hessian(ALag, [transpose(x) transpose(lambda)]); % hessian
    % Solve using pure newton's method
    [~, l_newton, l_iter_newton, l_err_newton] = pure_newton(ALag, grad_aLag, hess_aLag, [transpose(x) transpose(lambda)], l0, eps);
    results_L_newton = [results_L_newton; [l_newton,l_iter_newton, l_err_newton(end)]];
end

%% Print out table with results
mu = transpose(10.^(-4:10));
xStarNewton = results_L_newton(:, 1:N+nEqConst);
iterNewton = results_L_newton(:, N+nEqConst+1);
errorNewton = results_L_newton(:, N+nEqConst+2);
T_L_Newton = table(mu, xStarNewton, iterNewton, errorNewton)

%% Solve using SR1
results_L_SR1 = [];
for i = -4:2:10
    mu = 10^i;
    ALag = augmented_lagrangian(lambda, objFun, eqConst, mu); % define the augmented lagrangian for the current mu
    grad_aLag = gradient(ALag, [transpose(x) transpose(lambda)]); % gradient
    hess_aLag = hessian(ALag, [transpose(x) transpose(lambda)]); % hessian
    for gamma = 0:5:10
        H0_l = double(subs(hess_aLag,[transpose(x) transpose(lambda)],l0)) + gamma*eye(N+nEqConst);
        [L_x_SR1,l_x_SR1,l_iter_SR1,l_err_SR1] = SR1 (grad_aLag, H0_l, [transpose(x) transpose(lambda)], l0, eps);
        results_L_SR1 = [results_L_SR1; [l_x_SR1, l_iter_SR1,l_err_SR1(end)]];
    end
end

%% Print out table with results
mu = [];
gamma = [];
for i = -4:2:10
    for g = 0:5:10
        mu = [mu; 10^i];
        gamma = [gamma ; g];
    end
end
xStarSR1 = results_L_SR1(:, 1:N+nEqConst);
iterSR1 = results_L_SR1(:, N+nEqConst+1);
errorSR1 = results_L_SR1(:, N+nEqConst+2);
T_L_SR1 = table(mu, gamma, xStarSR1, iterSR1, errorSR1)