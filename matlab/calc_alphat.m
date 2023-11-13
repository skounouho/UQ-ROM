function alphat = calc_alphat(U1, delta, T)
    if size(U1, 2) == 1 % vector
        alphat = zeros(size(U1,1), T+1);
        alphat(:,1) = U1;
        alpha0 = U1; 
        I = eye(size(U1, 2));
        II = [I; zeros(size(U1, 2))];
        for i = 2:T+1
            AA = alpha0'*delta; S0 = delta'*delta;
            A = [alpha0, delta]; B=[AA, -S0; I, AA];
            t = (i-1)/T;
            alphat(:,i) = A * expm(t*B) * II * expm(-t*AA);
        end
    
    else % matrix
        alphat = zeros(size(U1,1), size(U1, 2), T+1);
        alphat(:,:,1) = U1;
        alpha0 = U1; 
        I = eye(size(U1, 2));
        II = [I; zeros(size(U1, 2))];
        for i = 2:T+1
            AA = alpha0'*delta; S0 = delta'*delta;
            A = [alpha0, delta]; B=[AA, -S0; I, AA];
            t = (i-1)/T;
            alphat(:,:,i) = A * expm(t*B) * II * expm(-t*AA);
        end
    end
end