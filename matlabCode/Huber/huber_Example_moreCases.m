% Huber Function - Example
% Code is adapted from Huber fitting code by Boyd https://web.stanford.edu/~boyd/papers/admm/
% Authors: Cassidy Buhler and Hande Benson
% This code generates random problems that are solved with Conjugate
% Gradient Method (CGM) with and without Cubic Regularization.
% We also include Boyd's ADMM for comparison.

% Select which model to run
model = 4; % 1 = CG with Powell Restarts, 2 = CG with Hybrid Cubic Reg, 3 = ADMM, 4= RUN ALL MODELS
names = ["CG WITH POWELL RESTARTS","CG WITH HYBRID CUBIC REGULARIZATION","ADMM"];

% Generate problem data
rng('default'); %setting seed
ms = [5000,10000,100000,100000];
ns = [2000, 2000,2000,5000];
numCases = length(ms);
numProbs = 100; %num problems

alpha = 1.0; % over-relaxation parameter (typical values between 1.0 and 1.8).
rho = 1.0; %augmented Lagrangian parameter

for i = 1:numCases
    m = ms(i);        % number of examples
    n = ns(i);       % number of features

    time = zeros(numProbs,3);
    iters = zeros(numProbs,3);
    status = zeros(numProbs,3);
    inCubic = zeros(numProbs,2);


    fprintf("Case %d\n", i);
    fprintf("m: %d, n: %d\n", m, n);

    for rr=1:numProbs

        fprintf("Case %d, Problem %d\n", i, rr);
        x0 = randn(n,1);
        A = randn(m,n);
        A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
        b = A*x0 + sqrt(0.01)*randn(m,1);
        b = b + 10*sprand(m,1,200/m);      % add sparse, large noise

        % Solve problem
        switch model
            case 1
                [x1, history1] = huber_cg_powellRestarts(A, b, alpha);
            case 2
                [x2, history2] = huber_cg_hybridCubic(A, b, alpha);
            case 3
                [x3, history3] = huber_admm(A, b, rho, alpha); %Boyd's ADMM code
                time(rr,:) =  history3.time; %saving time and iter per problem for performance profiles
                iters(rr,:) = history3.iters;
                status(rr,:) = history3.status;
            case 4
                fprintf([repmat('-',1,60),'\n']);
                fprintf("Case: %d, Problem %d\n",i, rr);
                fprintf("(%d,%d)\n", m, n);

                fprintf("\n%s\n",join(['--',names(1),'--'],''))
                [x1, history1] = huber_cg_powellRestarts(A, b, alpha);

                fprintf("\n%s\n",join(['--',names(2),'--'],''))
                [x2, history2] = huber_cg_hybridCubic(A, b, alpha);
                
                fprintf("\n%s\n",join(['--',names(3),'--'],''))
                [x3, history3] = huber_admm(A, b, rho, alpha); %Boyd's ADMM code

                 %saving time and iter per problem
                time(rr,:) = [history1.time, history2.time, history3.time];
                iters(rr,:) = [history1.iters, history2.iters, history3.iters];
                status(rr,:) = [history1.status, history2.status, history3.status];
                inCubic(rr,:) = [history1.powellRestart, history2.invokedCubic];
        end
    end  
end
