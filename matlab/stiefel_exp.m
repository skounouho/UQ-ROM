function U = stiefel_exp(U0, delta)
% U0: base point on St(N,n)
% delta: tangent vector on TuSt(N,n)

A = U0' * delta;
r = size(U0,2);
% thin qr-decomposition
[Q,R] = qr(delta - U0*A, 0);
% eigenvalue decomposition
[V,D] = eig([A, -R'; R, zeros(r)]);
MN = V * expm(D) * V' * [eye(r); zeros(r)];
M = real(MN(1:r, :));
N = real(MN(r+1:end, :));
U = U0*M + Q*N;

end