
N_iter = 10;
J = 1;
h = 1;
r = 4;

S_x = [0 1; 1 0];
S_z = [1 0; 0 -1];

C = 0;
hvals = h;
Jvals = J;
Cvals = C;
EperNvals = [];
for iter = 1:N_iter
  H = sparse(2^r,2^r);
  for pos = 1:r
    H = H + h*kron(speye(2^(pos-1)),kron(S_z,speye(2^(r-pos))));
    if pos ~= r
      H = H - J*kron(speye(2^(pos-1)),kron(S_x,kron(S_x,speye(2^(r-pos-1)))));
    end
  end
  [eigvec,eigval] = eigs(H,2,'smallestreal');
  eigval = diag(eigval);
  h = (eigval(2)-eigval(1))/2;
  EperNvals = [EperNvals,(eigval(1)+r*C)/r^iter];
  Sxmx1 = eigvec'*kron(speye(2^(r-1)),S_x)*eigvec;
  Sxmx2 = eigvec'*kron(S_x,speye(2^(r-1)))*eigvec;
  if any(abs(Sxmx1(:)-Sxmx2(:)) > 1e-12)
    error('')
  end
  J = J*Sxmx1(1,2)^2;
  hvals = [hvals,h];
  Jvals = [Jvals,J];
  C = r*C + (eigval(2)+eigval(1))/2;
  Cvals = [Cvals,C];
end

