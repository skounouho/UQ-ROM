function samples = dirichlet_sample(N_samples, Vs, K, gamma, U0)
% N_samples: number of samples to be generated
% X: initial dataset, in (d, N). they are tangential vectors on tangential
% space in particular
% K: number of samples in each iteration for dirichlet distribution
% gamma: scale in exponential distance
    d = size(Vs, 1);
    N = size(Vs, 2);
    samples = zeros(d, N_samples);
    for i = 1:N_samples
        y = Vs(:, randsample(1:N, K));
        y0 = y(:, 1);
        alpha = zeros(1, K);
        for j = 1:K
            yj = y(:,j);
            dist = calc_dist_metric_cano(yj, y0, U0);
            alpha(j) = exp(-gamma*dist*dist);
        end
%         disp(alpha);
        w = gamrnd(alpha,1);
        w = w ./ repmat(sum(w,2),1,K);
        y_interp = w*y';
        samples(:,i) = y_interp;
    end

end

   
    