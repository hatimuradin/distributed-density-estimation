function [log_pseudo_likelihood, logZ] = compute_log_pseudo_likelihood(SS, X, theta, exps, l_bound, r_bound)
    n = size(X,1);
    d = size(X,2);
    logZ = zeros(n,d);
    for i=1:n
        for j=1:d
            x_tmp = X(i,:);
            x_tmp(j) = 1;
            SS_nei = compute_SS( x_tmp, [], [], exps);
            buff = zeros(1,r_bound-l_bound+1);
            for v = l_bound:r_bound
                buff(v-l_bound+1) = (v.^(exps(j,:)) .* SS_nei') * theta;
            end
            logZ(i,j) = log_sum_exp(buff,2);
        end
    end
    log_pseudo_likelihood = d*theta'*SS - sum(sum(logZ));
end