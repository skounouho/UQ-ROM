function Delta2 = multi_parallel_translate(alphat, Delta)
    if size(Delta, 2) == 1
        T = size(alphat, 2) - 1;
        Delta2 = Delta;
        for i = T:1
            U2 = alphat(:, i);
            Delta2 = parallel_translate(U2, Delta2);
        end
    else
        T = size(alphat, 3) - 1;
        Delta2 = Delta;
        for i = fliplr(1:T)
            U2 = alphat(:, :, i);
            Delta2 = parallel_translate(U2, Delta2);
        end
        return
    end

end