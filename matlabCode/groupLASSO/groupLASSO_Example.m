% Group LASSO - Example
% Code is adapted from Group LASSO code by Boyd https://web.stanford.edu/~boyd/papers/admm/
% Authors: Cassidy Buhler and Hande Benson
% This code generates random problems that are solved with Conjugate
% Gradient Method (CGM) with and without Cubic Regularization.
% We also include Boyd's ADMM for comparison.

% Select which model to run
model = 4; % 1 = CG without Cubic Reg, 2 = CG with cubic Reg, 3 = ADMM, 4= RUN ALL MODELS
perfProf = true; % outputs a performance profile comparing all methods.
names = ["CG WITHOUT CUBIC REGULARIZATION","CG WITH CUBIC REGULARIZATION","ADMM"];

% Generate problem data
rng('default') %set seed
m = 100000;       % amount of data
N = 4;        % number of blocks
numProbs = 100; %num problems
K = 2000; % partition upperbound. randomly selects number between 1 and this number

if model ~= 4
    fprintf([repmat('-',1,60),'\n']);
    fprintf("%s\n", names(model)); % prints out solver name before solving problem(s)
    fprintf([repmat('-',1,60),'\n']);
end
alpha = 1.0; % over-relaxation parameter (typical values between 1.0 and 1.8).
rho = 1.0; %augmented Lagrangian parameter

time = zeros(numProbs,3);
iters = zeros(numProbs,3);
status = zeros(numProbs,3);
inCubic = zeros(numProbs,2);
for rr = 1:numProbs
    fprintf("Problem %d\n", rr);
    partition = randi(K, [N 1]);
    
    n = sum(partition); % number of features
    p = 100/n;          % sparsity density
    
    % generate block sparse solution vector
    x = zeros(n,1);
    start_ind = 1;
    cum_part = cumsum(partition);
    for i = 1:N
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
    
    lambdas = zeros(1,N);
    % lambda max
    start_ind = 1;
    for i = 1:N
        sel = start_ind:cum_part(i);
        lambdas(i) = norm(A(:,sel)'*b);
        start_ind = cum_part(i) + 1;
    end
    lambda_max = max(lambdas);
    % regularization parameter
    lambda = 0.01*lambda_max;
    
    xtrue = x;   % save solution
    % Solve problem
    switch model
        case 1
            [x1, history1] = groupLASSO_cg_noCubic(A, b, lambda, partition, alpha);
        case 2
            [x2, history2] = groupLASSO_cg_withCubic(A, b, lambda, partition, alpha);
        case 3
            [x3, history3] = groupLASSO_admm(A, b, lambda, partition, rho, alpha); %Boyd's ADMM code
        case 4 %run all models
            fprintf([repmat('-',1,60),'\n']);
            fprintf("%s\n",join(['--',names(1),'--'],''))% prints out solver name before solving problem
            [x1, history1] = groupLASSO_cg_noCubic(A, b, lambda, partition, alpha);
            fprintf("%s\n",join(['--',names(2),'--'],''))
            [x2, history2] = groupLASSO_cg_withCubic(A, b, lambda, partition, alpha);
            fprintf("%s\n",join(['--',names(3),'--'],''))
            [x3, history3] = groupLASSO_admm(A, b, lambda, partition, rho, alpha); %Boyd's ADMM code
            time(rr,:) = [history1.time, history2.time, history3.time]; %saving time and iterations per problem for performance profiles
            iters(rr,:) = [history1.iters, history2.iters, history3.iters];
            status(rr,:) = [history1.status, history2.status, history3.status];
            inCubic(rr,:) = [history1.powellRestart, history2.invokedCubic];
    end
    fprintf([repmat('-',1,60),'\n']);
end

if model == 4 && perfProf == true
    getPerformanceProfiles_invokeCubic(time,iters,status,inCubic) %performance profile of both CG's when cubic was invoked
    getPerformanceProfiles(time,iters,status) %performance profile for all 3
end
