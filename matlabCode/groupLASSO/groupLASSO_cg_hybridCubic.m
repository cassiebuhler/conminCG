function [z, history] = groupLASSO_cg_hybridCubic(A, b, lambda, p, alpha)

% Code is adapted from group LASSO code by Boyd https://web.stanford.edu/~boyd/papers/admm/
% Authors: Cassidy Buhler and Hande Benson

% [x, history] = groupLASSO_cg_hybridCubic(A, b, lambda, p, alpha)
%
% Solves the following problem via Conjugate Gradient WITH
% Cubic regularization:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda sum(norm(x_i))
%
% The input p is a K-element vector giving the block sizes n_i, so that x_i
% is in R^{n_i}.
%
% The solution is returned in the vector x.
%
% ﻿history is a struct that contains the objective values, l2 norm of gradients, time elapsed,
% number of iterations, solution status (0 = solved, 1 = Search direction is not descent direction,
% 2 = Iterations limit reached, 3 = search direction is undefined, 4 = Line search failed),
% and if cubic regularization was invoked (TRUE/FALSE)
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).


t_start = tic;
% Global constants and defaults
QUIET    = 1;
MAX_ITER = 1000;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

% Data preprocessing
[~, n] = size(A);


% check that sum(p) = total number of elements in x
if (sum(p) ~= n)
    error('invalid partition');
end

% cumulative partition
cum_part = cumsum(p);

% CG solver
x = 0.1*ones(n,1);
c = grad(A, b, lambda, x, cum_part);

pertcnt = 0;
nrst = n;
restart = false;

inPowell = false;
dprat = 0.0;
ddprat = 0.0;
lam = 0;
k = 0;

if ~QUIET
    fprintf("-----------------------------------------------------\n")
    fprintf("Iter |		Obj Value     	Residual	| P  \n")
    fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
end

while (k < MAX_ITER)
    k = k + 1;
    xTx = dot(x,x);
    cTc = dot(c,c);

    % Check for convergence
    if ( sqrt(cTc) <= sqrt(n)*ABSTOL + RELTOL*sqrt(xTx) ) % inftol*inftol*max(1.0, xTx) )
        status = 0; %solved
        break;
    end
    if ~QUIET
        f = objective(A, b, lambda, cum_part, x, x);
        fprintf("%4d :\t%14.6e\t %14.6e\t | %d \n", k, f, sqrt(cTc/max(1.0, xTx)),pertcnt);
    end
    % Compute step direction
    if ( restart == false )
        dx = -c;
    else
        % Test to see if the Powell restart criterion holds
        if ( nrst ~= n && restart > 1 && abs(dot(c,c0)/cTc) > 0.2 )
            if (pertcnt > 5) %powell restart
                dx = -c;
                inPowell = false;
                restart = false;
                pertcnt = 0;
            else

                if (nrst == 0)
                    nrst = n;
                end

                inPowell = true;
                if ( lam == 0.0 )
                    if ( dprat0 ~= 0.0 )
                        ddprat = ( abs(dot(c,c0)/cTc) - 2*prat + prat0 ) / ( (lam - olambda)*(oolambda-olambda) );
                    end
                    dprat0 = dprat;
                    if ( lam > 0 )
                        dprat = ( abs(dot(c,c0)/cTc) - prat ) / ( lam - olambda );
                    end
                    if (exist('prat','var'))
                        prat0 = prat;
                    end
                    prat = abs(dot(c,c0)/cTc);
                    alpha = alpha0;
                    restart = restart0;
                    x = x0;
                    dx = dx0;
                    c = c0;
                    c0 = c00;
                    if (exist('olambda','var'))
                        oolambda = olambda;
                    end
                    olambda = lam;
                    if ( lam == 0.0 )
                        lam = prat/0.2;
                    else
                        if ( pertcnt <= 2 )
                            lam = olambda - prat/dprat;
                        else
                            if ( dprat*dprat - 2*ddprat*prat > 0 && abs(ddprat) > 1e-8 )
                                lam = olambda + max( (-dprat + sqrt(dprat*dprat - 2*ddprat*prat))/ddprat, ...
                                    (-dprat - sqrt(dprat*dprat - 2*ddprat*prat))/ddprat );
                            else
                                lam = 2*olambda;
                            end
                        end
                        if ( lam < 1e-12 )
                            lam = 2*olambda;
                        end
                    end
                else
                    lam = 2*lam;
                end
                pertcnt = pertcnt + 1;
                k = k - 1;
            end
        else
            dx0 = dx;
            if ( pertcnt > 0 )
                lam = lam/2.0;
            else
                lam = 0.0;
            end
            dprat = 0.0;
            dprat0 = 0.0;
            prat0 = 0.0;
            pertcnt = 0;
        end

        % If performing a restart, update the Beale restart vectors
        if ( nrst == n )
            pt = alpha*dx;
            yt = c - c0;
            ytTyt = dot(yt,yt);
            cTyt = dot(pt,yt);
            cTct = dot(pt,pt);
        end

        p = alpha*dx;
        y = c - c0;
        pTc = dot(pt,c);
        yTc = dot(yt,c);

        u1 = -pTc/ytTyt;
        u2 = 2*pTc/cTyt - yTc/ytTyt;
        u3 = cTyt/ytTyt;

        if ( pertcnt == 0 )
            dx = -u3*c - u1*yt - u2*pt;
        else
            bracket = lam*lam + 2*lam*ytTyt/cTyt + ytTyt/cTct;
            a = -cTyt/(lam*cTyt + ytTyt);
            b1 = (-lam*cTyt - 2*ytTyt)*ytTyt/(cTyt*cTct*bracket*(lam*cTyt + ytTyt));
            d = lam/(bracket*(lam*cTyt + ytTyt));
            e = ytTyt/(cTct*bracket*(lam*cTyt + ytTyt));

            bracket = lam*lam+lam/u3+lam*ytTyt/cTyt+cTyt/(u3*cTct);
            denom = u3*u3*cTct*lam*bracket + u3*cTct*bracket;
            d = u3*u3*lam*(cTct/cTyt) / denom;

            dx = a*c + b1*pTc*pt + d*yTc*yt + e*yTc*pt + e*pTc*yt;
        end

        if ( nrst ~= n )
            if ( pertcnt == 0 )
                u1 = -dot(y,pt)/ytTyt;
                u2 = -dot(y,yt)/ytTyt + 2*dot(y,pt)/cTyt;
                u3 = dot(p,y);
                temp = cTyt/ytTyt*y + u1*yt + u2*pt;
                u4 = dot(temp, y);

                u1 = -dot(p,c)/u3;
                u2 = (u4/u3 + 1)*dot(p,c)/u3 - dot(c,temp)/u3;
                dx = dx - u1*temp - u2*p;
            else
                ptTy = dot(pt,y);
                ytTy = dot(yt,y);
                temp1 = -(a*y + b1*ptTy*pt + d*ytTy*yt + e*ytTy*pt + e*ptTy*yt);

                a2 = 1.0/(lam*u3+1.0);
                b2 = lam*b1;
                d2 = lam*d;
                e2 = lam*e;
                ptTp = dot(pt,p);
                ytTp = dot(yt,p);
                temp2 = a2*p + b2*ptTp*pt + d2*ytTp*yt + e2*ytTp*pt + e2*ptTp*yt;

                u10 = dot(temp1, c);
                u11 = dot(temp2, c);
                u12 = dot(temp1, y);
                u13 = dot(temp2, y);
                u14 = dot(temp2, p);
                u15 = dot(p, y);

                denom = lam*u14*(u15 + u12) + u13*u13;
                if denom < eps
                    status = 3;
                    fprintf("Search direction is undefined.\n")
                    break;
                end
                dx = dx + lambda*u14*u10*temp1/denom - ...
                    (u15 + u12)*u11*temp2/denom + ...
                    u13*u10*temp2/denom + u13*u11*temp1/denom;
            end
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
    if (exist('c0','var'))
        c00 = c0;
    end
    c0 = c;
    alpha0 = alpha;
    restart0 = restart;

    if ( restart == 0 )
        restart = 1;
    else
        if ( pertcnt == 0 )
            if ( nrst == n )
                nrst = 0;
            end
            nrst = nrst + 1;
            restart = 2;
        end
    end


    afind = @(a) objective(A, b, lambda, cum_part,  x + a*dx, x);
    [alpha,~,exitflag] = fminbnd(afind, 0, 10); % chose upperbound to balance speed and optimality

    if (exitflag ~= 1)
        status = 4;
        fprintf("Line search failed.\n");
        break;
    end

    % Take the step and update function value and gradient
    x = x0 + alpha*dx;
    c = grad(A, b, lambda, x, cum_part);
    history.objval(k)  = objective(A, b, lambda, cum_part, x, x);
    history.normGrad(k)  = norm(c);
end
if k == MAX_ITER
    status = 2;
    fprintf('Interations limit reached.\n')
end

elapsedTime = toc(t_start);
fprintf('Elapsed time is %f seconds.\n', elapsedTime);
fprintf('Iters = %d, invokedCubicReg = %s\n', k, string(inPowell == 1));

z = x;

history.time = elapsedTime;
history.iters = k;
history.status = status;
history.invokedCubic = inPowell==1;
end

function p = objective(A, b, lambda, cum_part, x, z)
obj = 0;
start_ind = 1;
for i = 1:length(cum_part)
    sel = start_ind:cum_part(i);
    obj = obj + norm(z(sel));
    start_ind = cum_part(i) + 1;
end
p = ( 1/2*sum((A*x - b).^2) + lambda*obj );
end


function c = grad(A, b, lambda, x, cum_part)
start_ind = 1;
c =A'*(A*x-b);
for i = 1:length(cum_part)
    sel = start_ind:cum_part(i);
    c(sel) = c(sel) +lambda*x(sel)/(norm(x(sel)));
    start_ind = cum_part(i) + 1;
end
end