function [f] = compute_momentpdfs(x, theta, exps, logZ)
    tmp = prod(repmat(x(:),1,size(exps,2)) .^ exps,1);
    P = exp(tmp * theta - logZ);
    f = tmp' .* P;
end