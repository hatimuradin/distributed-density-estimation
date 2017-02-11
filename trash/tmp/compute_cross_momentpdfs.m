function [f] = compute_cross_momentpdfs(x, theta, exps, logZ)
    tmp = prod(repmat(x(:),1,size(exps,2)) .^ exps,1);
    tmp2 = tmp' * tmp;
    P = exp(tmp * theta - logZ);
    f = tmp2(:) .* P;
end