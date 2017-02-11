% function X = generate_data(exps, theta, n)

% exps = [...
%     1 1 0 0 0 0 0 0 0 1 0 0;
%     1 0 1 1 1 0 0 0 0 1 1 0;
%     0 1 1 0 0 0 1 1 0 1 0 1;
%     0 0 0 1 0 1 1 0 1 0 1 1;
%     0 0 0 0 1 1 0 0 0 0 1 0;
%     0 0 0 0 0 0 0 1 1 0 0 1];

exps = [...
    1 1 0 0 0 0 0 0 0 1 0 0 0;
    1 0 1 1 1 0 0 0 0 1 1 0 1;
    0 1 1 0 0 0 1 1 0 1 0 1 1;
    0 0 0 1 0 1 1 0 1 0 1 1 1;
    0 0 0 0 1 1 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 1 1 0 0 1 0];

n = 1000;

% theta = 10^(-12) * [1;1;1;-1;1;1;1;1;1;-1;-1;-1];
theta = 10^(-12) * [1;1;1;-1;1;1;1;1;1;-1;-1;-1;-1];


    d = size(exps,1);
    tmp = ones(1,d);
    u_bound = 6;
    log_p = zeros(1,u_bound^d);
    
    tmp_X = zeros(u_bound^d, d);
    for k=1:u_bound^d
        tmp
        tmp_X(k,:) = tmp;
        log_p(k) = compute_poly(tmp, theta, exps);
        if k < u_bound^d
            j=find(tmp < u_bound, 1);
            tmp(1:j-1) = 1;
            tmp(j) = tmp(j) + 1;
        end
    end
    
    logZ = log_sum_exp(log_p,2);
    log_p = log_p - logZ;
    p = exp(log_p);
    cum_p = cumsum(p);
    
    plot(1:length(p),p);
    
    f = zeros(1,n);
    X = zeros(n,d);
    for i=1:n
        ff = find(rand() < cum_p, 1);
        if isempty(ff)
            ff = 0;
        end
        f(i) = ff+1;
        X(i,:) = tmp_X(f(i),:);
    end
    
    
    save X_fake_param_test.mat X;
% end