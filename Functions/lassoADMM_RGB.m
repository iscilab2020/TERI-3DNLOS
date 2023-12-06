function [zout, hist_r,hist_g,hist_b ] = lassoADMM_RGB(A, br,bg,bb, lambda, rho, alpha, MAX_ITER)
% lasso  Solve lasso problem via ADMM
%
% [z, history] = lasso(A, b, lambda, rho, alpha);
%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda || x ||_1
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%

if nargin==5
    MAX_ITER = 1000;
end
t_start = tic;
QUIET    = 1;

%Data preprocessing
[m, n] = size(A);

%ADMM solver
% x = zeros(n,1);
z = zeros(n,1);
u = zeros(n,1);

% cache the factorization
[L, U] = factor(A, rho);
% save a matrix-vector multiply
Atb = A'*[br,bg,bb];

if numel(lambda)==3
    lambda_r = lambda(1);
    lambda_g = lambda(2);
    lambda_b = lambda(3);
else
    lambda_r = lambda(1);
    lambda_g = lambda(1);
    lambda_b = lambda(1);
end
toc(t_start)
[zr, hist_r] = lassoADMMinternal(A, br, Atb(:,1), L, U, z, u, m, n, lambda_r, rho, alpha, MAX_ITER);
toc(t_start)
[zg, hist_g] = lassoADMMinternal(A, bg, Atb(:,2), L, U, z, u, m, n, lambda_g, rho, alpha, MAX_ITER);
toc(t_start)
[zb, hist_b] = lassoADMMinternal(A, bb, Atb(:,3), L, U, z, u, m, n, lambda_b, rho, alpha, MAX_ITER);
toc(t_start)

zout.r = zr;
zout.g = zg;
zout.b = zb;


if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end



if ~QUIET
    toc(t_start);
end
end

function [z, history] = lassoADMMinternal(A, b, Atb, L, U, z, u, m, n, lambda, rho, alpha, MAX_ITER)

%Global constants and defaults
QUIET    = 1;
% MAX_ITER = 1000;
ABSTOL   = 1e-8;
RELTOL   = 1e-8;


for k = 1:MAX_ITER

    % x-update
    q = Atb + rho*(z - u);    % temporary value
    if( m >= n )    % if skinny
       x = U \ (L \ q);
    else            % if fat
       x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
    end

    % z-update with relaxation
    zold = z;
    x_hat = alpha*x + (1 - alpha)*zold;
    z = shrinkage(x_hat + u, lambda/rho);

    % u-update
    u = u + (x_hat - z);

    % diagnostics, reporting, termination checks
    history.objval(k)  = objective(A, b, lambda, x, z);

    history.r_norm(k)  = norm(x - z);
    history.s_norm(k)  = norm(-rho*(z - zold));

    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);

    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end

end

end

function p = objective(A, b, lambda, x, z)
    p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) );
end

function z = shrinkage(x, kappa)
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end

function [L, U] = factor(A, rho)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A + rho*speye(n), 'lower' );
    else            % if fat
       L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end