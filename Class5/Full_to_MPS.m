
% --- get the ground state of the TFIM model for some L ---

TFIM_ED_OBC_sparse;     %this will provide us Psi, d, L

% --- Parameters ---
cut = L/2;    %




% --- transform Psi to an L-leg tensor

C = reshape(Psi,d*ones(1,L));
% Matlab specific issue:  kron() uses textbook Kronecker product that is
%                         compatible with 'Row-major' matrix layout.
%                         However, Matlab uses 'Column-major' storage:
%                         if we used kron() to build up our Hilbert space
%                         and operators, we need to reverse the legs after
%                         reshaping... 
%                         In python, there is no need for such
%                         transformation, as python uses row-major storage.
C = permute(C,L:-1:1);

A = {};
B = {};
% --- generate the left-matrices
if cut ~= 0
  C = reshape(C,[d,d^(L-1)]);
  [U,Lambda,V] = svd(C,'econ');    % econ makes Lambda a square matrix, i.e.,
                                   % trivially zero singular values are removed
  A{end+1} = U;
  C = reshape(Lambda*V',[size(Lambda,1),d*ones(1,L-1)]);
  for left = 2:cut
    sizeC1 = size(C,1); 
    sizeC2 = size(C,2);
    C = reshape(C,[sizeC1*sizeC2,d^(L-left)]);
    [U,Lambda,V] = svd(C,'econ');
    A{end+1} = reshape(U,[sizeC1,sizeC2,size(Lambda,1)]);
    if left ~= L
      C = reshape(Lambda*V',[size(Lambda,1),d*ones(1,L-left)]);
    else
      C = reshape(Lambda*V',[size(Lambda,1),1]);
    end
  end
  if cut == L 
    CoreTensor = Lambda;   %this will be actually 1 if Psi is normalized
  end
  CoreDim = size(C,1);
else
  CoreDim = 1;   
end
% --- generate the right matrices
if cut ~= L
  C = reshape(C,[CoreDim*d^(L-cut-1),d]);
  [U,Lambda,V] = svd(C,'econ');
  B{end+1} = V';
  C = reshape(U*Lambda,[CoreDim,d*ones(1,L-cut-1),size(Lambda,2)]);
  for right=2:(L-cut)
    sizeClast = size(C,L-cut+3-right);
    sizeCbeforelast=size(C,L-cut+2-right);
    C = reshape(C,[CoreDim*d^(L-cut-right),sizeCbeforelast*sizeClast]);
    [U,Lambda,V] = svd(C,'econ');
    B{end+1} = reshape(V',[size(Lambda,2),sizeCbeforelast,sizeClast]);
    if right ~= L
      C = reshape(U*Lambda,[CoreDim,d*ones(1,L-cut-right),size(Lambda,2)]);
    else
      C = reshape(U*Lambda,[CoreDim,size(Lambda,2)]);
    end
  end
  CoreTensor = Lambda;
  if cut ~= 0
    A{end} = tensorprod(A{end},U,3,1);
  end
end

% Check backwards  (we contract the MPS and compare with Psi;
for left = 1:cut
  if left == 1
    tmp = A{1};
  else
    tmp = tensorprod(tmp,A{left},left,1);
  end
end
if cut > 0
  tmp = tensorprod(tmp,CoreTensor,cut+1,1);
else
  tmp = CoreTensor;
end
for right = (L-cut):-1:1
  if cut > 0 || right < L-cut
    tmp = tensorprod(tmp,B{right},L-right+1,1);
  else
    tmp = tensorprod(tmp,B{right},L-right+1,1,"NumDimensionsA",1);   % we need to handle fake one dimensional leg in case of fully right canonical matrix
  end
end
Psi_test = reshape(permute(tmp,L:-1:1),[d^L,1]);
disp(['The difference norm is ',num2str(norm(Psi-Psi_test))]);


  



