
N_iter = 10;
J = 1;
h = 4;

S_x = [0 1; 1 0];
S_z = [1 0; 0 -1];

C = 0;
hvals = h;
Jvals = J;
Cvals = C;
for iter = 1:N_iter
  H = h*(kron(S_z,eye(2))) + h*kron(eye(2),S_z) - J*kron(S_x,S_x);
  [eigvec,eigval] = eig(H);
  eigval = diag(eigval);
  eigvec = eigvec(:,1:2);
  h = (eigval(2)-eigval(1))/2;
  Sxmx1 = eigvec'*(kron(S_x,eye(2)))*eigvec;
  Sxmx2 = eigvec'*(kron(eye(2),S_x))*eigvec;
  if any(abs(Sxmx1(:)-Sxmx2(:)) > 1e-12)
    error('')
  end
  J = J*Sxmx1(1,2)^2;
  hvals = [hvals,h];
  Jvals = [Jvals,J];
  C = 2*C + (eigval(2)+eigval(1))/2;
  Cvals = [Cvals,C];
end

