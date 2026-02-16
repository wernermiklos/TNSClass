% --- Simple Demo of the Lanczos Algorithm
%        We use the XXZ Model (PBC) Hamiltonian
% 
% --- Model Parameters ----
L = 20;
J = 1;
Delta = 0.5;
N_ITER = 200;  %number of Lanczos iterations.

% ---- 1-site operators ----
Sx_1site = sparse(0.5*[0 1;...
                       1 0]);
Sy_1site = sparse(0.5*[0 -1i;...
                       1i 0]);
Sz_1site = sparse(0.5*[1 0;...
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
for i = 1:L
  j = mod(i,L)+1;
  H = H + J*(Sx{i}*Sx{j} + Sy{i}*Sy{j} + Delta*Sz{i}*Sz{j});
end

%============================%
% ---- Lanczos iteration ----%
%============================%

PSI = zeros(2^L,N_ITER);
alpha = zeros(1,N_ITER);
beta = zeros(1,N_ITER);
PSI(:,1) = randn(2^L,1); 
PSI(:,1) = PSI(:,1) / norm(PSI(:,1));
wprime = H*PSI(:,1);
alpha(1) = wprime'*PSI(:,1);
w = wprime - alpha(1)*PSI(:,1);

H_Block_Lanczos = zeros(N_ITER,N_ITER);
H_Block_Lanczos(1,1) = alpha(1);
for j = 2:N_ITER
  disp(j);
  beta(j) = norm(w);
  if beta(j) == 0
    error('This should practically never happen.')
  end
  PSI(:,j) = w/beta(j);
  wprime = H*PSI(:,j) - beta(j)*PSI(:,j-1);
  alpha(j) = wprime'*PSI(:,j);
  w = wprime - alpha(j)*PSI(:,j);
  H_Block_Lanczos(j,j) = alpha(j);
  H_Block_Lanczos(j-1,j) = beta(j);
  H_Block_Lanczos(j,j-1) = beta(j);
end

ScalarProd_Check = PSI'*PSI;
H_Block_True = PSI'*(H*PSI);




