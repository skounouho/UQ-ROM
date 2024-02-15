function [frechet_mean,errs] = calc_frechet_mean_mat(samples, U0, eps)
    err = inf;
    c = 1;
    [d, r, N_samples] = size(samples);
    count = 0;
    errs = [];
    while err > eps
        V_mean = zeros(d, r);
        for i=1:N_samples
            V_mean = V_mean + real(stiefel_log(U0, samples(:,:,i), 1e-3));
        end
        V_mean = c * V_mean / N_samples;
        U0 = stiefel_exp(U0, V_mean);
        % calculate error
        err = norm(V_mean, 'fro');
        errs = [errs, err];
        count = count + 1;
        disp(['count = ', num2str(count), ', error = ', num2str(err)]);
%         pause();
    end

    frechet_mean = U0;
end