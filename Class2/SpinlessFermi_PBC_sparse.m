% --- Model Parameters ----
L = 20;
t = 1;
V = 0.5;

% ---- 1-site operators ----
cdag_1site = sparse([0 0;...
                     1 0]);
c_1site = cdag_1site';   %dagger

n_1site = cdag_1site*c_1site;

Id_1site = speye(2);

ph_1site = sparse([1 0; 0 -1]);

% ---- L-site spin operators ----
cdag = cell(1,L);
c = cell(1,L);
n = cell(1,L);
Ntot = [0];
ph_left = [1];

for pos = 1:L
  cdag{pos} = kron(ph_left,kron(cdag_1site,speye(2^(L-pos))));
  c{pos} = kron(ph_left,kron(c_1site,speye(2^(L-pos))));
  n{pos} = kron(speye(2^(pos-1)),kron(n_1site,speye(2^(L-pos))));
  Ntot = kron(Ntot,speye(2)) + kron(speye(2^(pos-1)),n_1site);    %iterative construction
  ph_left = kron(ph_left,ph_1site);
end

% ---- Hamiltonian ----
H = sparse(2^L,2^L);
for i = 1:L
  j = mod(i,L)+1;
  H = H - t*(cdag{i}*c{j} + cdag{j}*c{i}) + V*n{i}*n{j};
end

% ---- Diagonalization for the smallest 6 eigenvalues ----

[V,D] = eigs(H,6,'smallestreal');   %we ask for the 6 smallest eigenvalues

% -- PROBLEM:  particle number?

