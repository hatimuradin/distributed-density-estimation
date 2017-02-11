function [logP, logZ] = compute_log_likelihood(SS, n, theta, exps, l_bound, r_bound)

        d = size(exps,1);

        log_func = @(x) compute_poly(x,theta, exps);
        logZ = log_disc_integral_exp(log_func, l_bound, r_bound, d);
        
        logP = SS' * theta - n*logZ;
end