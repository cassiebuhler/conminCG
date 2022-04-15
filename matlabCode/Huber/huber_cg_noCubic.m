function [z, history] = huber_cg_noCubic(A, b, alpha)

% Code is adapted from Huber fitting code by Boyd https://web.stanford.edu/~boyd/papers/admm/
% Authors: Cassidy Buhler and Hande Benson

% This code minimizes the Huber loss function using Conjugate Gradient
% WITHOUT Cubic regularization.
%
% [z, history] =  huber_cg_withCubic(A, b, alpha);
%
% Solves the following problem via CG with cubic regularization:
%
%   minimize h( Az-b ) where h is the Huber loss function
% The solution is returned in the vector z.
%
% ﻿history is a struct that contains the objective values, l2 norm of gradients, time elapsed,
% number of iterations, solution status (0 = solved, 1 = Search direction is not descent direction, 
% 2 = Iterations limit reached, 3 = search direction is undefined, 4 = Line search failed), 
% and if a Powell restart was needed (TRUE/FALSE)
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).


t_start = tic;

% Global constants and defaults
QUIET    = 0;
MAX_ITER = 1000;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;
% Data preprocessing
[m, n] = size(A);


% CG solver
x = 10*ones(n,1);
c = grad(A, b, x, m, 1.0);

nrst = n;
restart = false;
inPowell = false;

for k = 1:MAX_ITER
    
    xTx = dot(x,x);
    cTc = dot(c,c);
    
    % Check for convergence
    if ( sqrt(cTc) <= sqrt(n)*ABSTOL + RELTOL*sqrt(xTx) ) 
        status = 0;
        break;
    end
    
    
    % Compute step direction
    if ( restart == false )
        dx = -c;
    else
        % Test to see if the Powell restart criterion holds
        if ( (nrst ~= n) && (abs(dot(c, c0)/cTc) > 0.2) )
            nrst = n;
            inPowell = true;
        end
        
        % If performing a restart, update the Beale restart vectors
        if ( nrst == n )
            pt = alpha*dx;
            yt = c - c0;
            ytTyt = dot(yt,yt);
            cTyt = dot(pt,yt);
        end
        
        p = alpha*dx;
        y = c - c0;
        pTc = dot(pt,c);
        yTc = dot(yt,c);
        
        u1 = -pTc/ytTyt;
        u2 = 2*pTc/cTyt - yTc/ytTyt;
        u3 = cTyt/ytTyt;
        dx = -u3*c - u1*yt - u2*pt;
        
        if ( nrst ~= n )
            u1 = -dot(y,pt)/ytTyt;
            u2 = -dot(y,yt)/ytTyt + 2*dot(y,pt)/cTyt;
            u3 = dot(p,y);
            temp = cTyt/ytTyt*y + u1*yt + u2*pt;
            u4 = dot(temp, y);
            
            u1 = -dot(p,c)/u3;
            u2 = (u4/u3 + 1)*dot(p,c)/u3 - dot(c,temp)/u3;
            dx = dx - u1*temp - u2*p;
        end
    end
    
    % Check that the search direction is a descent direction
    dxTc = dot(dx, c);
    if ( dxTc > 0 )
        status = 1;
        fprintf("CUBIT: Search direction is not a descent direction.\n");
        break;
    end
    
    % Save the current point
    x0 = x;
    c0 = c;
    
    if ( restart == 0 )
        restart = 1;
    else
        if ( nrst == n )
            nrst = 0;
        end
        nrst = nrst + 1;
        restart = 2;
    end
    
    afind = @(a) objective(A, b, x + a*dx);
    [alpha,~,exitflag] = fminbnd(afind, 0, 10);
    
    if (exitflag ~= 1)
        status = 4;
        fprintf("Line search failed.\n");
        break;
    end
    
    % Take the step and update function value and gradient
    x = x0 + alpha*dx;
    c = grad(A, b, x, m, 1.0);
    history.objval(k)  = objective(A, b, x);
    history.normGrad(k)  = norm(c);  
end

if k == MAX_ITER
    status = 2;
    fprintf('Interations limit reached.\n')
end
if ~QUIET
    elapsedTime = toc(t_start);
    fprintf('Elapsed time is %f seconds.\n', elapsedTime);
    fprintf('Iters = %d, PowellRestartNeeded = %s\n', k, string(inPowell == 1));
    
end
z = x;

history.time = elapsedTime;
history.iters = k;
history.status = status;
history.powellRestart = inPowell==1;
end

function p = objective(A, b, x)
p =  1/2*huber(A*x - b, 1.0);
end

function d = huber(x, gamma)
x1 = -x - 0.5;
x2 = 0.5*x.*x;
x3 = x - 0.5;
d = sum(x1.*(x <= -gamma) + x2.*(-gamma < x).*(x < gamma ) + x3.*(x >= gamma));
end

function c = huber_grad(x, n, gamma)
c = -ones(n,1).*(x <= -gamma) + x.*( -gamma < x).*(x < gamma ) + ones(n,1).*(x >= gamma);
end

function c = grad(A, b, x, n, gamma)
c = A'*huber_grad((A*x-b),n,gamma);
end
