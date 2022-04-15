% Huber Function - Example
% Code is adapted from Huber fitting code by Boyd https://web.stanford.edu/~boyd/papers/admm/
% Authors: Cassidy Buhler and Hande Benson
% This code generates random problems that are solved with Conjugate
% Gradient Method (CGM) with and without Cubic Regularization.
% We also include Boyd's ADMM for comparison.

diary huber_output.txt 

% Select which model to run
model = 4; % 1 = CG without Cubic Reg, 2 = CG with cubic Reg, 3 = ADMM, 4 = ALL MODELS
perfProf = true; % outputs a performance profile comparing all methods.
names = ["CG WITHOUT CUBIC REGULARIZATION","CG WITH CUBIC REGULARIZATION","ADMM"];

% Generate problem data
rng('default'); %setting seed
ms = [5000,10000,100000,100000];
ns = [2000, 2000,2000,5000];
numCases = length(ms);
numProbs = 100; %num problems

time = zeros(numProbs*numCases,3);
iters = zeros(numProbs*numCases,3);
status = zeros(numProbs*numCases,3);
inCubic = zeros(numProbs*numCases,2);

for i = 1:numCases
    m = ms(i);        % number of examples
    n = ns(i);       % number of features
    
    
    if model ~= 4
        fprintf([repmat('-',1,40),'\n']);
        fprintf("%s\n", names(model)); % prints out solver name before solving problem(s)
        fprintf([repmat('-',1,40),'\n']);
    end
    alpha = 1.0; % over-relaxation parameter (typical values between 1.0 and 1.8).
    rho = 1.0; %augmented Lagrangian parameter
    
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
                time(numProbs*(i-1)+rr,:) = [history1.time, history2.time, history3.time]; %saving time and iter per problem for performance profiles
                iters(numProbs*(i-1)+rr,:) = [history1.iters, history2.iters, history3.iters];
                status(numProbs*(i-1)+rr,:) = [history1.status, history2.status, history3.status];
                inCubic(numProbs*(i-1)+rr,:) = [history1.powellRestart, history2.invokedCubic];
        end
        fprintf([repmat('-',1,60),'\n']);
        
    end
end
if model == 4 && perfProf == true
    getPerformanceProfiles_invokeCubic(time,iters,status,inCubic) %performance profile of both CG's when cubic was invoked 
    getPerformanceProfiles(time,iters,status) %performance profile for all 3 
    
end
% save results 

% writetable(array2table(time,'VariableNames',{'NoCubic','WithCubic','ADMM'}),'huber_time.csv')
% writetable(array2table(iters,'VariableNames',{'NoCubic','WithCubic','ADMM'}),'huber_iters.csv')
% writetable(array2table(status,'VariableNames',{'NoCubic','WithCubic','ADMM'}),'huber_status.csv')
% writetable(array2table(inCubic,'VariableNames',{'NoCubic','WithCubic'}),'huber_inCubic.csv')

diary off 
