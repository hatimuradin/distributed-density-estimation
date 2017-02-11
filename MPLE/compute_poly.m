function [P] = compute_poly(X, theta, exps)
    P = prod(repmat(X(:),1,size(exps,2)) .^ exps,1) * theta;
end

