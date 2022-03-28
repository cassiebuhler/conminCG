% Huber Function - Example
% Code is adapted from Huber fitting code by Boyd https://web.stanford.edu/~boyd/papers/admm/
% Authors: Cassidy Buhler and Hande Benson
% Date last modified: March 24, 2022
% This code generates random problems that are solved with Conjugate
% Gradient Method (CGM) with and without Cubic Regularization.
% We also include Boyd's ADMM for comparison.
clear all; close all; clc;
% Select which model to run
model = 4; % 1 = CG without Cubic Reg, 2 = CG with cubic Reg, 3 = ADMM, 4 = ALL MODELS
names = ["CG WITHOUT CUBIC REGULARIZATION","CG WITH CUBIC REGULARIZATION","ADMM"];

% Generate problem data
rng('default'); %setting seed
m = 1000;        % number of examples
n = 500;       % number of features
numProbs = 100; %num problems

if model ~= 4
    fprintf([repmat('-',1,40),'\n']);
    fprintf("%s\n", names(model));
    fprintf([repmat('-',1,40),'\n']);
end
alpha = 1.0; % over-relaxation parameter (typical values between 1.0 and 1.8).
rho = 1.0; %augmented Lagrangian parameter

time = zeros(numProbs,3);
iters = zeros(numProbs,3);
ratios_time = zeros(numProbs,3);
ratios_iters = zeros(numProbs,3);

for rr=1:numProbs
    fprintf("Problem %d\n", rr);
    x0 = randn(n,1);
    A = randn(m,n);
    A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
    b = A*x0 + sqrt(0.01)*randn(m,1);
    b = b + 10*sprand(m,1,200/m);      % add sparse, large noise
    
    % Solve problem
    switch model
        case 1
            [x1, history1] = huber_cg_noCubic(A, b, alpha);
        case 2
            [x2, history2] = huber_cg_withCubic(A, b, alpha);
        case 3
            [x3, history3] = huber_admm(A, b, rho, alpha); %Boyd's ADMM code
        case 4
            fprintf([repmat('-',1,60),'\n']);
            fprintf("%s\n",join(['--',names(1),'--'],''))
            [x1, history1] = huber_cg_noCubic(A, b, alpha);
            fprintf("%s\n",join(['--',names(2),'--'],''))
            [x2, history2] = huber_cg_withCubic(A, b, alpha);
            fprintf("%s\n",join(['--',names(3),'--'],''))
            [x3, history3] = huber_admm(A, b, rho, alpha); %Boyd's ADMM code
            time(rr,:) = [history1.time, history2.time, history3.time];
            iters(rr,:) = [history1.iters, history2.iters, history3.iters];
            ratios_time(rr,:) = getRatio(time(rr,:),0); %compute ratio for performance profile
            ratios_iters(rr,:) = getRatio(iters(rr,:),1);
    end
    fprintf([repmat('-',1,60),'\n']);


end

if model == 4 % if ran all models, plot the performance prof
    ppIters_cgNoCubic = performanceProf(ratios_iters(:,1));
    ppIters_cgCubic = performanceProf(ratios_iters(:,2));
    ppIters_admm = performanceProf(ratios_iters(:,3));

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
    plot(ppIters_admm(:,1),ppIters_admm(:,2),'linewidth',2)
    legend({'CG without Cubic Reg','CG with Cubic Reg','ADMM'},'location','southeast','fontsize',14)
    ylim([0,1])
    
end

function ratios = getRatio(vec, isIter)
best = min(vec);
ratios = vec./best; %compute ratios
if isIter
    ratios(vec == 1000) = 1e8; %if it doesn't solve, set ratio to large value
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
