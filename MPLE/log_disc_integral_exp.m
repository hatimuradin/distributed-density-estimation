function log_intgrl = log_disc_integral_exp(log_func, l_bound, r_bound, d)
    disp('log_disc_integral_exp  start');
    L = (r_bound - l_bound + 1)^d;
    x = l_bound * ones(d,1);
    
    buff_size = 10000000;
    buff = zeros(buff_size,1);
    log_intgrl = -Inf;
    
    i = 0;
    for t=1:L
        if mod(t,10000) == 0
            fprintf('.');
        end
        for j=1:d
            if x(j) < r_bound
                x(j) = x(j)+1;
                x(1:j-1) = l_bound;
                break;
            end
        end
        i = i+1;
        
        buff(i) = log_func(x);
        if i == buff_size
            log_intgrl = log_sum_exp([log_intgrl; buff],1);
            i = 0;
        end
    end
    
    if i > 0
        log_intgrl = log_sum_exp([log_intgrl; buff(1:i)],1);
    end
    
    fprintf('\n');
    
    disp('log_disc_integral_exp  finish');
end