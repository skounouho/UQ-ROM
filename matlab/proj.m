function delta = proj(U, Z)
    % subroutine 2 in Bryner 2017
    I = eye(size(U,1));
    delta = U*(U'*Z-Z'*U)/2+(I - U*U')*Z;
end