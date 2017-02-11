function [theta, log_likelihood, logZ, iter] = fit_logpoly(SS, n, exps, theta_init, l_bound, r_bound)

%     alg = 'gradient';
    alg = 'Newton';

    M = length(SS);
    d = size(exps,1);
    
    wb = waitbar(0,'Optimizing...');
    theta = theta_init;
    [current_log_likelihood, current_logZ] = compute_log_likelihood(SS, n, theta, exps, l_bound, r_bound);

    MaxIter = 250;
    for iter = 1:MaxIter
        %%
        waitbar(iter/MaxIter,wb,['Optimizing... Iteration ' num2str(iter) ' of ' num2str(MaxIter)]);
        %% ESS

        
        func = @(x) compute_momentpdfs(x,theta, exps, current_logZ);
        ESS = disc_integral(func,l_bound, r_bound, d, M);
        
    %% Gradient Ascent
       grad = (SS - n*ESS);
       delta_x = grad / n;

    %% Newton method
    %     % Hessian
        if isequal(alg,'Newton')
            func = @(x) compute_cross_momentpdfs(x,theta, exps, current_logZ);
            H = disc_integral(func,l_bound, r_bound, d, M*M);
            H = reshape(H,M,M);
            H = -n*(H - ESS * ESS');
%             H = vpa(H);
            H = mp(H);
            grad = (SS - n*ESS);
            
            disp('inverting start');
            delta_x = (-grad'/H')';
            disp('inverting finish');

            disp('test precision start');
%             tmp = sum(sum(abs(inv(H) * H - eye(size(H,1)))));
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
              if lambda2/2 < (1e-6) % * n
                  disp(['converged in ' num2str(iter) ' steps.']);
                  break;
              end
            end
%             if ~isempty(find(H*delta_x + grad, 1))
%                 disp(['OH! ' num2str(sum(H*delta_x + grad))]);
%             end
        end
        %% Backtrack Line search
        lambda = 1;
        alpha = 0.49; beta = 0.5;

        delta_x = double(delta_x);
        while true
            [next_log_likelihood, next_logZ] = compute_log_likelihood(SS, n, theta + lambda*delta_x, exps, l_bound, r_bound);
            if next_log_likelihood >= current_log_likelihood + alpha*lambda*grad'*delta_x
                break;
            end
            lambda = lambda * beta;
%             lambda
%             next_log_likelihood
%             current_log_likelihood
%             current_log_likelihood + alpha*lambda*grad'*delta_x
        end
        
        if (next_log_likelihood <= current_log_likelihood)
            disp(['converged in ' num2str(iter) ' steps.']);
            disp('STOPPED BY THE SECOND STOPPING CRITERIA');
            break;
        end
        
        theta = theta + lambda * delta_x;

        current_log_likelihood = next_log_likelihood;
        current_logZ = next_logZ;
        
        current_log_likelihood
        
    end
    close(wb);
    
   log_likelihood = current_log_likelihood; 
   logZ = current_logZ;
end