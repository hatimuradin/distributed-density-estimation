clear all;
clc;
close all hidden;

addpath('AdvanpixMCT-3.9.9.11157');

% digits(50);
mp.Digits(50);

[X, X_v] = get_data();   %% get the train data (X) and validation data (X_v)
X = X(1:1000,:);

n = size(X,1);
d = size(X,2);

k = 8;        %% upper bound for sum of the power degrees of all dimension in each potential fuction
r = 3;        %% maximum size of the cluques

[SS, ~, exps] = compute_SS(X,k,r);  %% the sufficient statistic (SS) and a matrix (exps) that each column of it shows the power of different dimensions in one potential function.

M = length(SS);
% M

theta_init = zeros(M,1);

l_bound = min(unique(X)); %% with this assumption that all dimensions has the same set of allowed values
r_bound = max(unique(X)); %% with this assumption that all dimensions has the same set of allowed values
% l_bound = 1;
% r_bound = 6;

% [theta, log_likelihood, logZ, iter] = fit_logpoly(SS, n, exps, theta_init, l_bound, r_bound);
theta = fit_logpoly_mple(X,SS,exps,theta_init,l_bound,r_bound);

% [X, X_v] = get_data();
n = size(X,1);
SS = compute_SS(X,[],[],exps);

[log_likelihood, logZ] = compute_log_likelihood(SS, n, theta, exps, l_bound, r_bound);
vSS = compute_SS(X_v,k,r,exps);
vn = size(X_v,1);
[v_log_likelihood, v_logZ] = compute_log_likelihood(vSS,vn,theta,exps,l_bound,r_bound);
% save(['./results/X-5-6/' num2str(r) '-' num2str(k) '.mat'], ...
%     'log_likelihood','logZ','v_log_likelihood', 'v_logZ', 'iter');
save(['./results/X-5-6/' num2str(r) '-' num2str(k) '_mple_n=1000.mat'], ...
    'theta','log_likelihood','logZ','v_log_likelihood', 'v_logZ');

fprintf('log_likelihood= %.2f\n',(log_likelihood/size(X,1)));
fprintf('v_log_likelihood= %.2f\n',(v_log_likelihood/size(X_v,1)));
