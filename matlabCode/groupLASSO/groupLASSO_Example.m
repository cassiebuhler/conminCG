% Group LASSO - Example
% Code is adapted from Group LASSO code by Boyd https://web.stanford.edu/~boyd/papers/admm/
% Authors: Cassidy Buhler and Hande Benson
% Date last modified: March 11, 2022
% This code generates random problems that are solved with Conjugate
% Gradient Method (CGM) with and without Cubic Regularization.
% We also include Boyd's ADMM for comparison.

% Select which model to run
model = 4; % 1 = CG without Cubic Reg, 2 = CG with cubic Reg, 3 = ADMM, 4= RUN ALL MODELS
names = ["CG WITHOUT CUBIC REGULARIZATION","CG WITH CUBIC REGULARIZATION","ADMM"];

% Generate problem data
rng('default') %set seed
m = 1000;       % amount of data
K = 4;        % number of blocks

if model ~= 4
    fprintf([repmat('-',1,60),'\n']);
    fprintf("%s\n", names(model));
    fprintf([repmat('-',1,60),'\n']);
end
alpha = 1.0; % over-relaxation parameter (typical values between 1.0 and 1.8).
rho = 1.0; %augmented Lagrangian parameter


for rr = 1:100
    fprintf("Problem %d\n", rr);
    partition = randi(1000, [K 1]);
    
    n = sum(partition); % number of features
    p = 100/n;          % sparsity density
    
    % generate block sparse solution vector
    x = zeros(n,1);
    start_ind = 1;
    cum_part = cumsum(partition);
    for i = 1:K
        x(start_ind:cum_part(i)) = 0;
        if( rand() < p)
            % fill nonzeros
            x(start_ind:cum_part(i)) = randn(partition(i),1);
        end
        start_ind = cum_part(i)+1;
    end
    
    % generate random data matrix
    A = randn(m,n);
    
    % normalize columns of A
    A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n);
    
    % generate measurement b with noise
    b = A*x + sqrt(0.001)*randn(m,1);
    
    % lambda max
    start_ind = 1;
    for i = 1:K
        sel = start_ind:cum_part(i);
        lambdas(i) = norm(A(:,sel)'*b);
        start_ind = cum_part(i) + 1;
    end
    lambda_max = max(lambdas);
    
    % regularization parameter
    lambda = 0.001*lambda_max;
    xtrue = x;   % save solution
    % Solve problem
    switch model
        case 1
            [x, history] = groupLASSO_cg_noCubic(A, b, lambda, partition, alpha);
        case 2
            [x, history] = groupLASSO_cg_withCubic(A, b, lambda, partition, alpha);
        case 3
            [x, history] = groupLASSO_admm(A, b, lambda, partition, rho, alpha); %Boyd's ADMM code
        case 4 %run all models
            fprintf([repmat('-',1,60),'\n']);
            fprintf("%s\n",join(['--',names(1),'--'],''))
            [x1, history1] = groupLASSO_cg_noCubic(A, b, lambda, partition, alpha);
            fprintf("%s\n",join(['--',names(2),'--'],''))
            [x2, history2] = groupLASSO_cg_withCubic(A, b, lambda, partition, alpha);
            fprintf("%s\n",join(['--',names(3),'--'],''))
            [x3, history3] = groupLASSO_admm(A, b, lambda, partition, rho, alpha); %Boyd's ADMM code

    end
    
   fprintf([repmat('-',1,60),'\n']);
 
end
