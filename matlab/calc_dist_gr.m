function y = calc_dist_gr(U, V)
    s = svd(V' * U);
    y = sqrt(sum(acos(s).^2));
    y = real(y);
end