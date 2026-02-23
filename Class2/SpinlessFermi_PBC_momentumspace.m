% --- Model Parameters ----
L = 18;
t = 1;
V = 2;
N_elec = L/2;   %electron number

% ---- 1-site operators ----
cdag_1mode = sparse([0 0;...
                     1 0]);
c_1mode = cdag_1mode';   %dagger

n_1mode = cdag_1mode*c_1mode;

Id_1mode = speye(2);

ph_1mode = sparse([1 0; 0 -1]);

% ---- L-site spin operators ----
cdag = cell(1,L);
c = cell(1,L);
n = cell(1,L);
Ntot = sparse([0]);
mtot_tmp = sparse([0]);
ph_left = [1];

%k = (2pi/L*(m-1))

for m = 1:L
  cdag{m} = kron(ph_left,kron(cdag_1mode,speye(2^(L-m))));
  c{m} = kron(ph_left,kron(c_1mode,speye(2^(L-m))));
  n{m} = kron(speye(2^(m-1)),kron(n_1mode,speye(2^(L-m))));
  Ntot = kron(Ntot,speye(2)) + kron(speye(2^(m-1)),n_1mode);    %iterative construction
  mtot_tmp = kron(mtot_tmp,speye(2)) + kron(speye(2^(m-1)),n_1mode*(m-1));   %total momentum (measured in 'm') 
                                                                             %it's tmp because modulo L was not taken yet
  ph_left = kron(ph_left,ph_1mode);
end

% ---- particle number projector ----
Ntot_values = full(diag(Ntot));
States_selected = find(Ntot_values == N_elec);   %states having N_elec particles
MySecDim = length(States_selected);
Proj_N = sparse(1:MySecDim,States_selected',ones(1,MySecDim),...
                MySecDim,2^L);

% ---- total momentum operator ----
mtot = mod(Proj_N*(mtot_tmp*Proj_N.'),L);
mtot_values = full(diag(mtot));




% ---- formation of op2_cdag_c
op2_cdag_c = cell(L,L);
for m2 = 1:L
  c_times_Proj_N_T = c{m2}*Proj_N.';
  for m2tild = 1:L
    op2_cdag_c{m2tild,m2} = Proj_N*(cdag{m2tild}*c_times_Proj_N_T);   
  end
end



% ---- Hamiltonian  (for N_elec particles)
H = sparse(MySecDim,MySecDim);



for m1 = 1:L
  disp(m1);
  H = H - 2*t*cos(2*pi/L*(m1-1))*op2_cdag_c{m1,m1};
  for m2 = 1:L
    for m2tild = 1:L
      m1tild = mod(m1+m2-m2tild-1,L)+1;   %quasimomentum-conservation
      H = H + V* exp(1i*2*pi/L*(m2-m2tild))/L *op2_cdag_c{m1tild,m1}*op2_cdag_c{m2tild,m2};
    end
  end
end

% ---- Diagonalization for the smallest 6 eigenvalues per momentum ----

spectr = zeros(6,L);
for m = 1:L
  disp(m);
  States_m =(mtot_values == m-1);
  Hproj = H(States_m,States_m);
  spectr(:,m) = eigs(Hproj,6,'smallestreal');
end

 plot(0:2*pi/L:2*pi-2*pi/L,spectr)


