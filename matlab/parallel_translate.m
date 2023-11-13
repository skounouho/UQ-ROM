function Delta2 =  parallel_translate(U2, Delta)
l = norm(Delta);
Delta1 = proj(U2, Delta);
Delta2 = Delta1 * l / norm(Delta1);
end