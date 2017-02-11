function [X, X_v, P, H] = get_data()

%%
%     n = 10000;
%     d = 5;
%     
%     X = [ mvnrnd( zeros(d,1), eye(d), n/2 ); mvnrnd( 4*ones(d,1), 2*eye(d), n/2 )];
%     
%     X = X - repmat(min(X,[],1),n,1);
%     X = 10 * X ./ repmat(max(X,[],1),n,1);
%     
%     X = floor(X);
%     X(X == 10) = 9;
%     
%     X = X + 1;
%     save('X-5-10.mat','X');

%%%%%%%%%%%%%%%%%%%%%%%%%

%     n_train = 10000;
%     n_validation = 10000;
%     n = n_train + n_validation;
%     d = 5;
%     
%     X = [ mvnrnd( zeros(d,1), eye(d), n/2 ); mvnrnd( 4*ones(d,1), 2*eye(d), n/2 )];
%     X = X(randperm(n),:);
%     
%     X = X - repmat(min(X,[],1),n,1);
%     X = 6 * X ./ repmat(max(X,[],1),n,1);
%     
%     X = floor(X);
%     X(X == 6) = 5;
%     
%     X = X + 1;
%     X_v = X(n_train+1:end,:);
%     X = X(1:n_train,:);
%     save('X-5-6.mat','X','X_v');
% 

%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     n_train = 10000;
%     n_validation = 10000;
%     n = n_train + n_validation;
%     d = 5;
% 
%     l_bound = 1;
%     r_bound = 6;
%     L = (r_bound - l_bound + 1)^d;
%     z = l_bound * ones(d,1);
%     P = zeros(L,1);
%     mu = [zeros(d,1) , 4*ones(d,1), [1;1;1;4;4], [4;4;1;1;1]];
% %     mu = [zeros(d,1) , 4*ones(d,1) , [5;3;2;4;2]];
%     Sigma(:,:,1) = eye(d);
%     Sigma(:,:,2) = 2*eye(d);
%     Sigma(:,:,3) = eye(d);
%     Sigma(:,:,4) = 2*eye(d);
% %     Sigma(:,:,3) = [...
% % 5 4 0 3 0
% % 4 5 1 2 0
% % 0 1 5 1 3
% % 3 2 1 5 1
% % 0 0 3 1 5];
%     
%     for t=1:L
%         if mod(t,10000) == 0
%             fprintf('.');
%         end
%         for j=1:d
%             if z(j) < r_bound
%                 z(j) = z(j)+1;
%                 z(1:j-1) = l_bound;
%                 break;
%             end
%         end
%         for q=1:size(mu,2)
%             P(t) = P(t) + mvnpdf(z,mu(:,q),Sigma(:,:,q))^2;
%         end
%     end
%     P = P / sum(P);
% 
%     h = mnrnd(n,P);
%     X = zeros(n,d);
%     z = l_bound * ones(d,1);
%     i = 1;
%     for t=1:L
%         if mod(t,10000) == 0
%             fprintf('.');
%         end
%         for j=1:d
%             if z(j) < r_bound
%                 z(j) = z(j)+1;
%                 z(1:j-1) = l_bound;
%                 break;
%             end
%         end
%         X(i:(i+h(t)-1),:) = repmat(z',h(t),1);
%         i = i + h(t);
%     end
%     
%     X = X(randperm(n),:);
%     X_v = X(n_train+1:end,:);
%     X = X(1:n_train,:);
%     
%     ind = P ~= 0;
%     H = P(ind)' * log(P(ind));
%     
%     save('X-5-6.mat','X','X_v','P','H');

    load('X-5-6.mat');
    
end