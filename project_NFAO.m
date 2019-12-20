% NFAO 
% Computational Assignment - Version 5
% Kilian & Sara

%% define format
format long
%% Initialize figures numbering
fig = 1;
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
eigenvalues = eig(hess_star) % eigenvalues are non negative
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

x0 = [0.9, 0.9, 1, 1.1]; % row vector of dimension N
eps = 1e-6; % epsilon
tic()
[X_newton, x_newton,iter_newton,err_newton] = pure_newton (f,grad_f,hess_f, x, x0, eps);
x_newton, iter_newton
grad_newton = double(subs(grad_f,x,x_newton))
elapsed_time_d = toc()

%% Question f :

tic()
x0 = [0.9, 0.9, 1, 1.1];
gamma = 0;
H0 = double(subs(hess_f,x,x0)) + gamma*eye(N);
[X_SR1,x_SR1,iter_SR1,err_SR1] = SR1 (grad_f,H0, x, x0, eps);
x_SR1, iter_SR1
grad_SR1 = double(subs(grad_f,x,x_SR1))
elapsed_time_f = toc()

%% Question g : log plot
figure(fig)
fig = fig+1;
xx = 1:iter_newton-1;
yy = err_newton(xx);
plot(xx,log(yy), 'x')
a = polyfit(xx,log(yy),1);
hold on
plot(xx,a(1)*xx+a(2), 'Color',[0.3 0.3 0.3])
title(['Newton log ||gradf(xk)|| slope = ',num2str(a(1))])
xlabel('iterations k')
ylabel('log ||gradf(x_k)||')

figure(fig)
fig = fig+1;
xx = 1:iter_SR1-1;
yy = err_SR1(xx);
plot(xx,log(yy),'x')
b = polyfit(xx,log(yy),1)
hold on
plot(xx,b(1)*xx+b(2), 'Color',[0.3 0.3 0.3])
title(['SR1 log ||gradf(xk)|| slope = ',num2str(b(1))])
xlabel('iterations k')
ylabel('log ||gradf(x_k)||')


%% Question g : averaging with n different x0
%rng(43);
n = 20;
variance = 0.08; %x0 close from x_star
fitn = [0,0];
fits = [0,0];
for k = 1:n
    x_star = [1,1,1,1];
    gauss_noise = variance*randn(1,4);
    x0 = x_star + gauss_noise;
    H0 = double(subs(hess_f,x,x0));
    [~,~,itN,errN] = pure_newton (f,grad_f,hess_f, x, x0, eps);
    [~,~,itS,errS] = SR1 (grad_f,H0, x, x0, eps);
    xn = 1:(itN-1);
    xs = 1:(itS-1);
    
    fitn = fitn + polyfit(xn,log(errN(xn)),1);
    fits = fits + polyfit(xs,log(errS(xs)),1);
end
R_newton = -fitn(1)/n
R_sr1 = -fits(1)/n


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
KKT_sol = solve(kkt_conditions, [transpose(x) transpose(lambda) transpose(mu)], 'Real', true)
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
% solve gradient of the lagrangian in order to [x; lambda] = 0
%% Solve using Newton's method
z0 = randi(10, size(x,1) + size(lambda,1), 1) % random initial guess (entries are integers between 0 and 10)
eps = 1e-6; % considered zero
[zk_newton, iter_newton, err_newton] = newton(gradient(lag, [x; lambda]), z0, [x; lambda], eps)
% zk_newton is the concatenation of x_star and lambda_star found by newton's method
%% Solve using Broyden's method (5 trials)
z0Broyden = [];
iterBroyden = [];
zStarBroyden = [];
for i = 1:5
    z0 = randi(10,size(x,1) + size(lambda,1), 1); %random initial guess (entries are integers between 1 and 10)
    z0Broyden = [z0Broyden; transpose(z0)];
    eps = 1e-8; %epsilon
    A0 = rand(size(x,1) + size(lambda,1)); % generate a random squared matrix with the same dimension of the kkt system
    A0 = A0*transpose(A0); % make it positive semidefinite
    [zk_broyden, iter_broyden, err_broyden] = broyden(gradient(lag, [x;lambda]), z0, A0, [x; lambda], eps);
    iterBroyden = [iterBroyden; iter_broyden];
    zStarBroyden = [zStarBroyden; transpose(zk_broyden)]
end
%% Print out table with results (Broyden's method)
T_Broyden = table(z0Broyden, iterBroyden, zStarBroyden)

%% Question e (Quadratic penalty method)
% definition of the quadratic penalty function in quadratic_penalty.m
x0 = rand(N, 1); % set a random initial point
mu0 = 1; % initial mu value
eps = 1e-6; % considered 0
%% Quadratic penalty method using SR1 to solve the minimization problems
[X_Q_SR1, x_Q_SR1, iter_Q_SR1, cond_HQ_SR1] = quadratic_penalty_method_SR1(x, x0, objFun, eqConst, mu0, eps)
%% Plot the error over the iterations
% Compute the error of each iteration
e_k = repmat(transpose(x_star),size(X_Q_SR1,1),1)-X_Q_SR1; % vector of (x_star - x_k)
e_Q_SR1 = vecnorm(transpose(e_k)); % compute the norm of each entry of e_k

figure(fig)
fig = fig+1;
xx = 1:iter_Q_SR1;
yy = e_Q_SR1(xx);
plot(xx,log(yy), 'x')
a = polyfit(xx,log(yy),1);
hold on
plot(xx,a(1)*xx+a(2), 'Color',[0.3 0.3 0.3])
title(['SR1 Quadratic | slope = ',num2str(a(1))])
xlabel('iterations k')
ylabel('log ||x_\ast-x_k||')

%% Quadratic penalty method using Newton's method to solve minimization problems
[X_Q_Newton,x_Q_Newton, iter_Q_Newton, cond_HQ_Newton] = quadratic_penalty_method_Newton(x, x0, objFun, eqConst, mu0, eps)
%% Plot the error over the iterations
% Compute the error of each iteration
e_k = repmat(transpose(x_star),size(X_Q_Newton,1),1)-X_Q_Newton; % vector of (x_star - x_k)
e_Q_Newton = vecnorm(transpose(e_k)); % compute the norm of each entry of e_k

figure(fig)
fig = fig+1;
xx = 1:iter_Q_Newton;
yy = e_Q_Newton(xx);
plot(xx,log(yy), 'x')
a = polyfit(xx,log(yy),1);
hold on
plot(xx,a(1)*xx+a(2), 'Color',[0.3 0.3 0.3])
title(['Newton Quadratic | slope = ',num2str(a(1))])
xlabel('iterations k')
ylabel('log ||x_\ast-x_k||')

%% Question e (Augmented lagrangian method)
x0 = rand(N, 1); % set a random initial point
mu0 = 1; % initial mu value
lambda0 = rand(nEqConst, 1); % set a random initial lambda
eps = 1e-6; % considered 0
%% Augmented lagrangian method using SR1 to solve the minimization problems
[X_AL_SR1,x_AL_SR1, lambda_AL_SR1, iter_AL_SR1, cond_HAL_SR1] = augmented_lagrangian_method_SR1(x, x0, objFun, eqConst, mu0,lambda0, eps)
%% Plot the error over the iterations
% Compute the error of each iteration
e_k = repmat(transpose(x_star),size(X_AL_SR1,1),1)-X_AL_SR1; % vector of (x_star - x_k)
e_AL_SR1 = vecnorm(transpose(e_k)); % compute the norm of each entry of e_k
figure(fig)
fig = fig+1;
xx = 1:iter_AL_SR1;
yy = e_AL_SR1(xx);
plot(xx,log(yy), 'x')
a = polyfit(xx,log(yy),1);
hold on
plot(xx,a(1)*xx+a(2), 'Color',[0.3 0.3 0.3])
title(['SR1 Augmented Lagrangian | slope = ',num2str(a(1))])
xlabel('iterations k')
ylabel('log ||x_\ast-x_k||')

%% Augmented lagrangian method using Newton's method to solve the minimization problems
[X_AL_Newton, x_AL_Newton,lambda_AL_Newton, iter_AL_Newton, cond_HAL_Newton] = augmented_lagrangian_method_Newton(x, x0, objFun, eqConst, mu0,lambda0, eps)
%% Plot the error over the iterations
% Compute the error of each iteration
e_k = repmat(transpose(x_star),size(X_AL_Newton,1),1)-X_AL_Newton; % vector of (x_star - x_k)
e_AL_Newton = vecnorm(transpose(e_k)); % compute the norm of each entry of e_k

figure(fig)
fig = fig+1;
xx = 1:iter_AL_Newton;
yy = e_AL_Newton(xx);
plot(xx,log(yy), 'x')
a = polyfit(xx,log(yy),1);
hold on
plot(xx,a(1)*xx+a(2), 'Color',[0.3 0.3 0.3])
title(['Newton Augmented Lagrangian | slope = ',num2str(a(1))])
xlabel('iterations k')
ylabel('log ||x_\ast-x_k||')