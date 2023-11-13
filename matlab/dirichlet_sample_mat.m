function [samples, weights] = dirichlet_sample_mat(N_samples, robs, U0, K, gamma)
% This code aims for Dirichlet sampling on Stiefel manifold.
% N_samples: number of generated samples
% robs: the ROB candidates
% U0: the base ROB to generate tangential space
% K: number of samples in each iteration for Dirichlet sampling
% gamma: scale in exponential distance

    [n, r, N0] = size(robs); % [total dof, rank, number of ROB candidates]
    samples = zeros(n, r, N_samples);
    weights = zeros(N_samples, K);
%     f_idxs = zeros(N_samples, K);

    % calculate the distante matrix
    rob_dist = zeros(N0);
    for i = 1:N0
        for j = 1:N0
            rob_dist(i,j) = calc_dist_metric_cano(robs(:,:,i), robs(:,:,j));
        end
    end

    % project robs to the tangential space
    Vs = zeros(size(robs)); % projected samples in tangent space
    tau = 1e-5;
    for i = 1:N0
        Ui = robs(:, :, i); 
        Vs(:,:,i) = real(stiefel_log(U0, Ui, tau));
    end

    % Dirichlet sampling
    for i = 1:N_samples
        idx = randsample(1:N0, K);
        y = Vs(:, :, idx);
        distances = rob_dist(idx, idx);
        alpha = zeros(1, K);
        for j = 1:K
            dist = distances(1, j);
            alpha(j) = exp(-gamma*dist*dist);
        end
        w = gamrnd(alpha,1); % gamma distribution
        w = w ./ repmat(sum(w,2),1,K); % beta distribution
        weights(i, :) = w;
        delta = zeros(n, r); % the new generated tangential vector
        for j = 1:K
            delta = delta + w(j) * y(:,:,j);
        end
        samples(:,:,i) = stiefel_exp(U0, delta); % project back to Stiefel manifold
        if mod(i, 500) == 0
            disp(['Generated ', num2str(i), ' samples']);
        end
    end

end