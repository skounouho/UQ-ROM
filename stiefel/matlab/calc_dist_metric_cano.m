function dist_cano = calc_dist_metric_cano(U, V, varargin)
    n = size(U,1);
    if ~isempty(varargin)
        % U and V are on tangential space, not Stiefel space
        W0 = varargin{1};
        dist_cano = sqrt(trace(U'*(eye(n)-1/2*(W0*W0'))*V));
    else
        % U and V are on Stiefel spaces
        delta = real(stiefel_log(U, V, 1e-4)); % project V to the tangential space of U
        dist_cano = sqrt(trace(delta'*(eye(n)-1/2*(U*U'))*delta));
    end
    dist_cano = real(dist_cano);
end