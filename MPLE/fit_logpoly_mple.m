function theta = fit_logpoly_mple(X,SS,exps,theta_init,l_bound,r_bound)

%     X = X(1,:);

    n = size(X,1);
    
    %     alg = 'sgd';
    alg = 'Newton';

    d = size(exps,1);
    M = size(exps,2);
    k = max(max(exps));
    
    wb = waitbar(0,'Optimizing...');
    
    
    theta = theta_init;
%     load('theta.mat');

    [current_log_pseudo_likelihood, current_logZ] = compute_log_pseudo_likelihood(SS, X, theta, exps, l_bound, r_bound);

    MaxIter = 250;
    for iter = 1:MaxIter
        tic
        waitbar(iter/MaxIter,wb,['Optimizing... Iteration ' num2str(iter) ' of ' num2str(MaxIter)]);

%         save theta.mat theta
        
    %% Newton method
        if isequal(alg,'Newton')    

            grad = d*SS;
            H = zeros(M,M);
            for i=1:n
%                 if mod(i,1000) == 0
%                     i
%                 end
                for j=1:d
                    x_tmp = X(i,:);
                    x_tmp(j) = 1;
                    SS_nei = compute_SS( x_tmp, [], [], exps);
                    mmnt = zeros(1,2*k+1);
                    for pw = 0:2*k
                       buff = zeros(1,r_bound-l_bound+1);
                       for v = l_bound:r_bound
                            buff(v-l_bound+1) = v^pw * exp( (v.^(exps(j,:)) .* SS_nei') * theta - current_logZ(i,j));
                       end
                       mmnt(pw+1) = sum(buff);
                    end
                    ESS = SS_nei .* mmnt(exps(j,:)+1)';
                    grad = grad - ESS;
                    tmp = repmat(exps(j,:),M,1);
                    tmp = tmp + tmp';
                    H = H - ( ((SS_nei * SS_nei') .* mmnt(tmp+1)) - ESS * ESS');
                end
            end
            toc
            tic
            
            H = mp(H);
%             save H.mat H
%             rank(H)
            
            disp('inverting start');
            delta_x = (-grad'/H')';
            disp('inverting finish');

            disp('test precision start');
            tmp = sum(abs(grad' + delta_x'*H));
            if tmp > 1e-10
                tmp
                error('Precision error');
            end
            disp('test precision finish');


            if isequal(alg,'Newton')
              lambda2 = grad'*delta_x;
              if lambda2 < 0
                  lambda2
                  error('lambda2 is less than zero!');
              end
              if lambda2/2 < (1e-2) * n
                  disp(['converged in ' num2str(iter) ' steps.']);
                  break;
              end
            end
%             if ~isempty(find(H*delta_x + grad, 1))
%                 disp(['OH! ' num2str(sum(H*delta_x + grad))]);
%             end
        end
        toc
        
        %% Backtrack Line search
        lambda = 1;
        alpha = 0.49; beta = 0.5;

        delta_x = double(delta_x);
        while true
            [next_log_pseudo_likelihood, next_logZ] = compute_log_pseudo_likelihood(SS, X, theta + lambda*delta_x, exps, l_bound, r_bound);
            if next_log_pseudo_likelihood >= current_log_pseudo_likelihood + alpha*lambda*grad'*delta_x
                break;
            end
            lambda = lambda * beta;
            lambda
        end
        
        if (next_log_pseudo_likelihood <= current_log_pseudo_likelihood)
            disp(['converged in ' num2str(iter) ' steps.']);
            disp('STOPPED BY THE SECOND STOPPING CRITERIA');
            break;
        end
        
        theta = theta + lambda * delta_x;

        current_log_pseudo_likelihood = next_log_pseudo_likelihood;
        current_logZ = next_logZ;
        
        current_log_pseudo_likelihood
        
    end
    close(wb);
    
end