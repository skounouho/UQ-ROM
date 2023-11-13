function U2 = stiefel_exp_euclidean(U1, Delta)
    I = eye(size(U1, 2));
    II = [I; zeros(size(U1, 2))];
    AA = U1'*Delta; S0 = Delta'*Delta;
    A = [U1, Delta]; B=[AA, -S0; I, AA];
    t = 1;
    U2 = A * expm(t*B) * II * expm(-t*AA);
end

