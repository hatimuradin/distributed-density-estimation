function intgrl = disc_integral(func, l_bound, r_bound, d, M)
    disp('disc_integral  start');
    L = (r_bound - l_bound + 1)^d;
    
    intgrl = zeros(M,1);
    
    x = l_bound * ones(d,1);
    for t=1:L
        if mod(t,1000) == 0
            fprintf('.');
        end
        for j=1:d
            if x(j) < r_bound
                x(j) = x(j)+1;
                x(1:j-1) = l_bound;
                break;
            end
        end
        intgrl = intgrl + func(x);
    end
    fprintf('\n');
    disp('disc_integral  finish');
end