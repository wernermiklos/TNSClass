% --- Model Parameters ----
L = 12;
J = 1;
Delta = 4;

% ---- 1-site operators ----
Sx_1site = 0.5*[0 1;...
                1 0];
Sy_1site = 0.5*[0 -1i;...
                1i 0];
Sz_1site = 0.5*[1 0;...
                0 -1];
Id_1site = eye(2);

% ---- L-site spin operators ----
Sx = cell(1,L);
Sy = cell(1,L);
Sz = cell(1,L);

for pos = 1:L
  Sx{pos} = kron(eye(2^(pos-1)),kron(Sx_1site,eye(2^(L-pos))));   % id * id * ...* id * Sx * id * ... * id 
  Sy{pos} = kron(eye(2^(pos-1)),kron(Sy_1site,eye(2^(L-pos))));
  Sz{pos} = kron(eye(2^(pos-1)),kron(Sz_1site,eye(2^(L-pos))));
end

% ---- Hamiltonian ----
H = zeros(2^L);
for i = 1:L
  j = mod(i,L)+1;
  H = H + J*(Sx{i}*Sx{j} + Sy{i}*Sy{j} + Delta*Sz{i}*Sz{j});
end

% ---- Full diagonalization ----

[V,D] = eig(H);
scatter(diag(D),ones(2^L,1),'.')



