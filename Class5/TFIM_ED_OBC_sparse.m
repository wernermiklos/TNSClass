% --- Here this script is just used to generate some meaningful input ---

% --- Model Parameters ----
L = 16;
d = 2;
J = 1;
h = 1.1;

% ---- 1-site operators ----
Sx_1site = sparse([0 1;...
                   1 0]);     %0.5 factor removed to match the BRG code,
Sy_1site = sparse([0 -1i;...
                   1i 0]);
Sz_1site = sparse([1 0;...
                   0 -1]);
Id_1site = speye(2);

% ---- L-site spin operators ----
Sx = cell(1,L);
Sy = cell(1,L);
Sz = cell(1,L);

for pos = 1:L
  Sx{pos} = kron(speye(2^(pos-1)),kron(Sx_1site,speye(2^(L-pos))));   % id * id * ...* id * Sx * id * ... * id 
  Sy{pos} = kron(speye(2^(pos-1)),kron(Sy_1site,speye(2^(L-pos))));
  Sz{pos} = kron(speye(2^(pos-1)),kron(Sz_1site,speye(2^(L-pos))));
end

% ---- Hamiltonian ----
H = sparse(2^L,2^L);
for i = 1:(L-1)
  j = i+1;
  H = H - J*Sx{i}*Sx{j} + h*Sz{i};
end
H = H + h/2*Sz{L};

% ---- Full diagonalization ----

[V,D] = eigs(H,2,'smallestreal');   %we ask for the 2 smallest eigenvalues

Psi = V(:,1);   %this is the ground state

