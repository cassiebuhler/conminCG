% Huber Function - Example
% Code is adapted from Huber fitting code by Boyd https://web.stanford.edu/~boyd/papers/admm/
% Authors: Cassidy Buhler and Hande Benson
% Date last modified: March 11, 2022
% This code generates random problems that are solved with Conjugate
% Gradient Method (CGM) with and without Cubic Regularization.
% We also include Boyd's ADMM for comparison.

% Select which model to run
model = 4; % 1 = CG without Cubic Reg, 2 = CG with cubic Reg, 3 = ADMM, 4 = ALL MODELS
names = ["CG WITHOUT CUBIC REGULARIZATION","CG WITH CUBIC REGULARIZATION","ADMM"];

% Generate problem data
rng('default'); %setting seed
m = 10000;        % number of examples
n = 5000;       % number of features
if model ~= 4
    fprintf([repmat('-',1,40),'\n']);
    fprintf("%s\n", names(model));
    fprintf([repmat('-',1,40),'\n']);
end

for rr=1:100 %run 100 instances
    fprintf("Problem %d\n", rr);
    x0 = randn(n,1);
    A = randn(m,n);
    A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
    b = A*x0 + sqrt(0.01)*randn(m,1);
    b = b + 10*sprand(m,1,200/m);      % add sparse, large noise
    lambda_max = norm( A'*b, 'inf' );
    lambda = 0.001*lambda_max;
    
    % Solve problem
    switch model
        case 1
            [x, history] = huber_cg_noCubic(A, b, 1.0);
        case 2
            [x, history] = huber_cg_withCubic(A, b, 1.0);
        case 3
            [x, history] = huber_admm(A, b, 1.0, 1.0); %Boyd's ADMM code
        case 4
            fprintf([repmat('-',1,60),'\n']);
            fprintf("%s\n",join(['--',names(1),'--'],''))
            [x1, history1] = huber_cg_noCubic(A, b, 1.0);
            fprintf("%s\n",join(['--',names(2),'--'],''))
            [x2, history2] = huber_cg_withCubic(A, b, 1.0);
            fprintf("%s\n",join(['--',names(3),'--'],''))
            [x3, history3] = huber_admm(A, b, 1.0, 1.0); %Boyd's ADMM code
            
    end
fprintf([repmat('-',1,60),'\n']);

end

