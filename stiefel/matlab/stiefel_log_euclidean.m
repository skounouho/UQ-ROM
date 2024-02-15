function Delta = stiefel_log_euclidean(U1, U2, tau)
% The projection algorithm to the tangential space based on Euclidean
% metric
% U1: base point on Stiefel manifold, U2: target point on Stiefel manifold
% tau: the convergence threshold, usually set as 1e-3~1e-5
delta = 1; % coefficient to update Delta
J = 1000; % max iteration step
T = 20; % number of segments in the geodesic
err = tau+1; % initial error
dim = size(U1, 2);

Delta = U2 - U1;  % initialization
Delta = proj(U1, Delta);

alphat = calc_alphat(U1, Delta, T);
j=0;
while err > tau && j<J
    j = j+1;
    if dim == 1
        alpha1 = alphat(:, end); 
    else
        alpha1 = alphat(:, :, end); 
    end
    w = alpha1 - U2;
    w = proj(alpha1, w);
    w = multi_parallel_translate(alphat, w);
    err = norm(w);
    Delta = Delta - delta * w;
    alphat = calc_alphat(U1, Delta, T);
end




