function [ SS, S, exps] = compute_SS( X, k, r, exps_given)
    d = size(X,2);
    
    if nargin == 4
        exps = exps_given;
        n = size(X,1);
        S = zeros(n,size(exps,2));
        for i=1:n
            S(i,:) = prod(repmat(X(i,:)' , 1, size(exps,2)) .^ exps, 1);
        end
        SS = sum(S,1)';
        return;
    elseif nargin ~= 3
        error('Wrong number of inputs!');
    end
    
    S = X;
    L = S;
    exps = eye(d);
    last_exps = exps;
%     L
    for i=2:k
        W = [];
        new_exps = [];
        for j=1:d
            for t=1:size(L,2)
                tmp = last_exps(:,t);
                if sum(tmp) == k
                    continue;
                end
                tmp(j) = tmp(j)+1;
                if sum(tmp > 0) > r
                    continue;
                end
                W = [ W , X(:,j).*L(:,t) ];
                new_exps = [new_exps, tmp];
            end
        end
        [~,ia] = unique(new_exps','rows');
        ia = sort(ia);
        new_exps = new_exps(:,ia);
        W = W(:,ia);
        
        S = [S,W];
        L = W;
        last_exps = new_exps;
        exps = [exps, new_exps];
%         L
    end
    
    SS = sum(S,1)';
end