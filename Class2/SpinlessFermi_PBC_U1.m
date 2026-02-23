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

% ---- L-site spin operators (%full, non-symmetric versions) ----
cdag_full = cell(1,L);
c_full = cell(1,L);
n_full = cell(1,L);
Ntot_full = [0];
ph_left_full = [1];

disp('Formation of full ops:')
for pos = 1:L
  disp(pos);
  cdag_full{pos} = kron(ph_left_full,kron(cdag_1site,speye(2^(L-pos))));
  c_full{pos} = kron(ph_left_full,kron(c_1site,speye(2^(L-pos))));
  n_full{pos} = kron(speye(2^(pos-1)),kron(n_1site,speye(2^(L-pos))));
  Ntot_full = kron(Ntot_full,speye(2)) + kron(speye(2^(pos-1)),n_1site);    %iterative construction
  ph_left_full = kron(ph_left_full,ph_1site);
end

% ---- Sectorization ----
Ntot_values = full(diag(Ntot_full));
Ntot_values_unique = unique(Ntot_values);    %the list of possible charges
SecDim = length(Ntot_values_unique);   %the number of sectors
States_of_Sec = cell(1,SecDim);
for mysec = 1:SecDim
  States_of_Sec{mysec} = (Ntot_values ==  Ntot_values_unique(mysec));
end
c = cell(1,L);
cdag = cell(1,L);
n = cell(1,L);
disp('Sectorization of ops:')
for pos = 1:L
  disp(pos);
  c{pos} = cell(SecDim,SecDim);
  cdag{pos} = cell(SecDim,SecDim);
  n{pos} = cell(SecDim,SecDim);
  for sec_i = 1:SecDim
    States_i = States_of_Sec{sec_i};
    for sec_f = 1:SecDim
      States_f =States_of_Sec{sec_f};
      cblock = c_full{pos}(States_f,States_i);
      if nnz(cblock) > 0    %number of nonzero elements > 0
        c{pos}{sec_f,sec_i} = cblock;
      end
      cdagblock = cdag_full{pos}(States_f,States_i);
      if nnz(cdagblock) > 0    %number of nonzero elements > 0
        cdag{pos}{sec_f,sec_i} = cdagblock;
      end
      nblock = n_full{pos}(States_f,States_i);
      if nnz(nblock) > 0    %number of nonzero elements > 0
        n{pos}{sec_f,sec_i} = nblock;
      end
    end
  end
end


% ---- Hamiltonian ----
H = cell(SecDim,SecDim);
disp ('formation of H:')
for i = 1:L
  disp(i);
  j = mod(i,L)+1;
  H = SecMatAdd(H, ...
                SecMatDot(cdag{i},c{j},-t),1,1);
  H = SecMatAdd(H, ...
                SecMatDot(cdag{j},c{i},-t),1,1);
  H = SecMatAdd(H,...
                SecMatDot(n{i},n{j},V),1,1);
end

% ---- Diagonalization for the smallest 6 eigenvalues ----
Sec_eigvals = cell(1,SecDim);
disp('Sector-wise diagonalization:')
for sec = 1:SecDim
  disp(sec)
  Sec_eigvals{sec} = eigs(H{sec,sec},6,'smallestreal');   %we ask for the 6 smallest eigenvalues
end


% ------------------------- INTERNAL FUNCTIONS --------------------------- %


function out = SecMatDot(A,B,fact)
  if size(A,2) ~= size(B,1)
    error('Contracted SecDim mismatch')
  end
  out = cell(size(A,1),size(B,2));
  for sec_f = 1:size(A,1)
    for sec_i = 1:size(B,2)
      for sec_c = 1:size(A,2)
        Ablock = A{sec_f,sec_c};
        Bblock = B{sec_c,sec_i};
        if ~isempty(Ablock) && ~isempty(Bblock)
          out_block = fact*Ablock*Bblock;
          if isempty(out{sec_f,sec_i})
            out{sec_f,sec_i} = out_block;
          else
            out{sec_f,sec_i} = out{sec_f,sec_i}+out_block;
          end
        end
      end
    end
  end
end

function out = SecMatAdd(A,B,fact_A, fact_B)
  if size(A,1) ~= size(B,1) || size(A,2) ~= size(B,2)
    error('SecDim mismatch')
  end
  out = cell(size(A,1),size(B,2));
  for sec_f = 1:size(A,1)
    for sec_i = 1:size(B,2)
      Ablock = A{sec_f,sec_i};
      Bblock = B{sec_f,sec_i};
      if ~isempty(Ablock) || ~isempty(Bblock)
        if isempty(Ablock)
          out_block = fact_B*Bblock;
        elseif isempty(Bblock)
          out_block = fact_A*Ablock;
        else
          out_block = fact_A*Ablock + fact_B*Bblock;
        end
        if isempty(out{sec_f,sec_i})
          out{sec_f,sec_i} = out_block;
        else
          out{sec_f,sec_i} = out{sec_f,sec_i}+out_block;
        end
      end
    end
  end
end
