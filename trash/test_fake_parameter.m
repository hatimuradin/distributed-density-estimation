clear all;
clc;
close all hidden;

addpath('AdvanpixMCT-3.9.9.11157');

% digits(50);
mp.Digits(50);

% exps = [...
%     1 1 0 0 0 0 0 0 0 1 0 0;
%     1 0 1 1 1 0 0 0 0 1 1 0;
%     0 1 1 0 0 0 1 1 0 1 0 1;
%     0 0 0 1 0 1 1 0 1 0 1 1;
%     0 0 0 0 1 1 0 0 0 0 1 0;
%     0 0 0 0 0 0 0 1 1 0 0 1];

% theta_true = [1;1;1;-1;1;1;1;1;1;-1;-1;-1];

%X = generate_data(exps,theta_true,1000);
load('X_fake_param_test.mat');

X = X(:,[2,3,4]);
-5.3735e+03
exps = [...
      1 1 0;
      1 0 1;
      0 1 1];

n = size(X,1);
d = size(X,2);


SS = compute_SS(X,[],[],exps);

M = length(SS);

theta_init = zeros(M,1);

l_bound = min(unique(X));
r_bound = max(unique(X));
% l_bound = 1;
% r_bound = 6;

[theta, log_likelihood, logZ, iter] = fit_logpoly(SS, n, exps, theta_init, l_bound, r_bound);
%theta = fit_logpoly_mple(X,SS,exps,theta_init,l_bound,r_bound);

%n = size(X,1);
%SS = compute_SS(X,[],[],exps);

%[log_likelihood, logZ] = compute_log_likelihood(SS, n, theta, exps, l_bound, r_bound);
%vSS = compute_SS(X_v,k,r,exps);
%vn = size(X_v,1);
%[v_log_likelihood, v_logZ] = compute_log_likelihood(vSS,vn,theta,exps,l_bound,r_bound);
% save(['./results/X-5-6/' num2str(r) '-' num2str(k) '.mat'], ...
%     'log_likelihood','logZ','v_log_likelihood', 'v_logZ', 'iter');
%save(['./results/X-5-6/' num2str(r) '-' num2str(k) '_mple_n=1000.mat'], ...
%     'theta','log_likelihood','logZ','v_log_likelihood', 'v_logZ');

%fprintf('log_likelihood= %.2f\n',(log_likelihood/size(X,1)));
%fprintf('v_log_likelihood= %.2f\n',(v_log_likelihood/size(X_v,1)));
