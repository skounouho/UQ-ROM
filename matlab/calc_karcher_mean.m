function [karcher_mean,errs] = calc_karcher_mean(samples, U0, eps, tau)
    err = inf;
    if length(size(samples)) == 2
        samples = reshape(samples, size(samples,1), 1, size(samples,2));
    end
    [d, r, N_samples] = size(samples);
    count = 0;
    errs = [];
    while err > eps
        V_mean = zeros(d, r);
        for i=1:N_samples
            V_mean = V_mean + stiefel_log_euclidean(U0, samples(:,:,i), tau);
        end
        V_mean = V_mean / N_samples;
        U0 = stiefel_exp_euclidean(U0, V_mean);
        % calculate error
        err = norm(V_mean, 'fro');
        errs = [errs, err];
        count = count + 1;
        disp(['count = ', num2str(count), ', error = ', num2str(err)]);
    end

    karcher_mean = U0;
end