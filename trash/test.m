load X.mat;

n = size(X,1);
d = size(X,2);

k = 6;
r = 3;

[SS, ~, exps] = compute_SS(X,k,r);

M = length(SS);
M

theta = zeros(M,1);

l_bound = 1;
r_bound = 10;


% log_func = @(x) compute_poly(x,theta, exps);
% logZ = log_disc_integral_exp(log_func, l_bound, r_bound, d);
logZ = 11.5129;

% func = @(x) compute_cross_momentpdfs(x,theta, exps, logZ);
% H = disc_integral(func,l_bound, r_bound, d, M*M);

%%

disp('disc_integral  start');
L = (r_bound - l_bound + 1)^d;

intgrl = zeros(M*M,1);

x = l_bound * ones(d,1);
for t=1:L
    if mod(t,100) == 0
        fprintf('.');
    end
    for j=1:d
        if x(j) < r_bound
            x(j) = x(j)+1;
            x(1:j-1) = l_bound;
            break;
        end
    end
%     
%     tmp = prod(repmat(x(:),1,size(exps,2)) .^ exps,1);
%     tmp2 = tmp' * tmp;
%     P = exp(tmp * theta - logZ);
%     f = tmp2(:) .* P;
    f = ones(M*M,1);
    
    intgrl = intgrl + f;
end
fprintf('\n');
disp('disc_integral  finish');

%%
