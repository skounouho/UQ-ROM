function dist_euc = calc_dist_metric_euc(U, V)
    delta = stiefel_log(U, V, 1e-4);
    dist_euc = sqrt(trace(delta'*delta));
end