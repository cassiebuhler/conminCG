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
m = 10000;       % amount of data
K = 5;        % number of blocks
numProbs = 100; %num problems

if model ~= 4
    fprintf([repmat('-',1,60),'\n']);
    fprintf("%s\n", names(model));
    fprintf([repmat('-',1,60),'\n']);
end
alpha = 1.0; % over-relaxation parameter (typical values between 1.0 and 1.8).
rho = 1.0; %augmented Lagrangian parameter

time = zeros(numProbs,3);
iters = zeros(numProbs,3);
ratios_time = zeros(numProbs,3);
ratios_iters = zeros(numProbs,2);

for rr = 1:numProbs
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
            [x1, history1] = groupLASSO_cg_noCubic(A, b, lambda, partition, alpha);
        case 2
            [x2, history2] = groupLASSO_cg_withCubic(A, b, lambda, partition, alpha);
        case 3
            [x3, history3] = groupLASSO_admm(A, b, lambda, partition, rho, alpha); %Boyd's ADMM code
        case 4 %run all models
            fprintf([repmat('-',1,60),'\n']);
            fprintf("%s\n",join(['--',names(1),'--'],''))
            [x1, history1] = groupLASSO_cg_noCubic(A, b, lambda, partition, alpha);
            fprintf("%s\n",join(['--',names(2),'--'],''))
            [x2, history2] = groupLASSO_cg_withCubic(A, b, lambda, partition, alpha);
            fprintf("%s\n",join(['--',names(3),'--'],''))
            [x3, history3] = groupLASSO_admm(A, b, lambda, partition, rho, alpha); %Boyd's ADMM code
            time(rr,:) = [history1.time, history2.time, history3.time];
            iters(rr,:) = [history1.iters, history2.iters, history3.iters]; %exclude admm from iterations pp 
            ratios_time(rr,:) = getRatio(time(rr,:),0); %compute ratio for performance profile
            ratios_iters(rr,:) = getRatio(iters(rr,1:2),1);%exclude admm from iterations pp 

    end
    
   fprintf([repmat('-',1,60),'\n']);
 
end

if model == 4 % if ran all models, plot the performance prof
    ppIters_cgNoCubic = performanceProf(ratios_iters(:,1));
    ppIters_cgCubic = performanceProf(ratios_iters(:,2));

    ppTime_cgNoCubic = performanceProf(ratios_time(:,1));
    ppTime_cgCubic = performanceProf(ratios_time(:,2));
    ppTime_admm = performanceProf(ratios_time(:,3));
    
    figure
    title("Performance Profile: Time",'fontsize',18)
    xlabel("Tau",'fontsize',16)
    ylabel("Probability",'fontsize',16)
    hold on
    plot(ppTime_cgNoCubic(:,1),ppTime_cgNoCubic(:,2),'linewidth',2)
    plot(ppTime_cgCubic(:,1),ppTime_cgCubic(:,2),'linewidth',2)
    plot(ppTime_admm(:,1),ppTime_admm(:,2),'linewidth',2)
    legend({'CG without Cubic Reg','CG with Cubic Reg','ADMM'},'location','southeast','fontsize',14)
    ylim([0,1])
    
    
    figure
    title("Performance Profile: Iterations",'fontsize',18)
    xlabel("Tau",'fontsize',16)
    ylabel("Probability",'fontsize',16)
    hold on
    plot(ppIters_cgNoCubic(:,1),ppIters_cgNoCubic(:,2),'linewidth',2)
    plot(ppIters_cgCubic(:,1),ppIters_cgCubic(:,2),'linewidth',2)
    legend({'CG without Cubic Reg','CG with Cubic Reg'},'location','southeast','fontsize',14)
    ylim([0,1])
    
end

function ratios = getRatio(vec, isIter)
best = min(vec);
ratios = vec./best; %compute ratios
if isIter
    ratios(vec == 1000) = NaN; %if it doesn't solve, set ratio to nan 
end
end

function performanceData = performanceProf(ratios) %compute performance prof
ratio = sort(ratios);
uniqueVals = unique(ratio);
counts=[];
for i=1:length(uniqueVals)
    counts(i) = sum(ratio == uniqueVals(i));
end
probs = [];
taus = [];
prob = 0;
for j=1:length(counts)
    prob = prob + counts(j)/length(ratio);
    probs(j) = prob;
    taus(j) = uniqueVals(j);
end
performanceData = [taus;probs]';
end
