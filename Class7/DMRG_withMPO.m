L0 = 4;   %initial chain length
L = 100;
J = 1;
Delta = 2;
d = 2;   %local dimension
M = 32;
HamMPO = MPO_XXZ_full(L0,J,Delta);
Nsweep = 2;

% INFINITE CHAIN (CHAIN-GROWTH) PART
% initialization of 1-site left and right block using the first and last MPO
% matrix:
LeftBlocks = {reshape(HamMPO{1},[d,d,size(HamMPO{1},4)])};   %remove 1-dim leg
RightBlocks = {reshape(HamMPO{L0},[d,d,size(HamMPO{L0},3)])};  %remove 1-dim leg
% Chain-growth iteration:
disp('INFINITE CHAIN (CHAIN-GROWTH) ALGORITHM')
disp('L_act En M TrErr')
for L_act = L0:2:L
  dimLeft = size(LeftBlocks{end},1);
  dimRight = size(RightBlocks{end},1);

  % form [[left]*]
  newLeftBlock = reshape(permute(tensorprod(LeftBlocks{end},HamMPO{L_act/2},3,3,"NumDimensionsA",3),...
                                 [1,3,2,4,5]), [dimLeft*d,dimLeft*d, size(HamMPO{L_act/2},4)]);
  % form [*[right]]
  newRightBlock = reshape(permute(tensorprod(HamMPO{L_act/2+1},RightBlocks{end},4,3,"NumDimensionsA",4),...
                                 [1,4,2,5,3]), [dimRight*d,dimRight*d, size(HamMPO{L_act/2+1},3)]);

  % Diagonalization
  dimSB = size(newLeftBlock,1)*size(newRightBlock,1);
  myHdotPsi = @(x) internal_HdotPsi(x,newLeftBlock,newRightBlock);
  [EigVec,En] = eigs(myHdotPsi,dimSB,[],1,'smallestreal');
  EigVec = reshape(EigVec,[dimLeft*d,dimRight*d]);
  
  % Schmidt decomposition and truncation
  [U,S,Vconj] = svd(EigVec,'econ');
  Truncerr = 0;
  if size(S,1) > M
    U = U(:,1:M);
    Vconj = Vconj(:,1:M);
    S = diag(S);
    Truncerr = S(M+1:end)'*S(M+1:end);
    S = S(1:M);
  end
  fprintf('%d %d %d %e\n',L_act,En,size(S,1),Truncerr)
  V = conj(Vconj);

  % Renormalization
  LeftBlocks{end+1} = permute(tensorprod(tensorprod(U',newLeftBlock,2,1),U,2,1),[1,3,2]);
  RightBlocks{end+1} = permute(tensorprod(tensorprod(V',newRightBlock,2,1),V,2,1),[1,3,2]);
  % Updating the MPO with the two new sites
  HamMPO = [HamMPO(1:L_act/2),HamMPO{L_act/2}, HamMPO{L_act/2+1}, HamMPO{L_act/2+1:end}];
end

disp('SWEEPING')
disp('left En M TrErr')
Energies = [En];
RightBlocks(end) = [];
RightBlocks(end) = [];
for sweep_iter = 1:Nsweep
  
  %center --> right end
  while length(RightBlocks) > 1
    left = length(LeftBlocks);
    dimLeft = size(LeftBlocks{end},1);
    dimRight = size(RightBlocks{end},1);
    % form [[left]*]
    newLeftBlock = reshape(permute(tensorprod(LeftBlocks{end},HamMPO{left+1},3,3,"NumDimensionsA",3),...
                                   [1,3,2,4,5]), [dimLeft*d,dimLeft*d, size(HamMPO{L_act/2},4)]);
    % form [*[right]]
    newRightBlock = reshape(permute(tensorprod(HamMPO{left+2},RightBlocks{end},4,3,"NumDimensionsA",4),...
                                   [1,4,2,5,3]), [dimRight*d,dimRight*d, size(HamMPO{L_act/2+1},3)]);
    % Diagonalization
    dimSB = size(newLeftBlock,1)*size(newRightBlock,1);
    myHdotPsi = @(x) internal_HdotPsi(x,newLeftBlock,newRightBlock);
    [EigVec,En] = eigs(myHdotPsi,dimSB,[],1,'smallestreal');
    Energies(end+1) = En;
    EigVec = reshape(EigVec,[dimLeft*d,dimRight*d]);  
    % Schmidt decomposition and truncation
    [U,S,Vconj] = svd(EigVec,'econ');
    Truncerr = 0;
    if size(S,1) > M
      U = U(:,1:M);
      Vconj = Vconj(:,1:M);
      S = diag(S);
      Truncerr = S(M+1:end)'*S(M+1:end);
      S = S(1:M);
    end
    fprintf('%d %d %d %e\n',left,En,size(S,1),Truncerr)
    V = conj(Vconj);
    % Renormalization
    LeftBlocks{end+1} = permute(tensorprod(tensorprod(U',newLeftBlock,2,1),U,2,1),[1,3,2]);
    RightBlocks(end) = [];
  end

  % left end <-- right end
  while length(LeftBlocks) > 1
    left = length(LeftBlocks);
    dimLeft = size(LeftBlocks{end},1);
    dimRight = size(RightBlocks{end},1);
    % form [[left]*]
    newLeftBlock = reshape(permute(tensorprod(LeftBlocks{end},HamMPO{left+1},3,3,"NumDimensionsA",3),...
                                   [1,3,2,4,5]), [dimLeft*d,dimLeft*d, size(HamMPO{L_act/2},4)]);
    % form [*[right]]
    newRightBlock = reshape(permute(tensorprod(HamMPO{left+2},RightBlocks{end},4,3,"NumDimensionsA",4),...
                                   [1,4,2,5,3]), [dimRight*d,dimRight*d, size(HamMPO{L_act/2+1},3)]);
    % Diagonalization
    dimSB = size(newLeftBlock,1)*size(newRightBlock,1);
    myHdotPsi = @(x) internal_HdotPsi(x,newLeftBlock,newRightBlock);
    [EigVec,En] = eigs(myHdotPsi,dimSB,[],1,'smallestreal');
    EigVec = reshape(EigVec,[dimLeft*d,dimRight*d]);  
    % Schmidt decomposition and truncation
    [U,S,Vconj] = svd(EigVec,'econ');
    Truncerr = 0;
    if size(S,1) > M
      U = U(:,1:M);
      Vconj = Vconj(:,1:M);
      S = diag(S);
      Truncerr = S(M+1:end)'*S(M+1:end);
      S = S(1:M);
    end
    fprintf('%d %d %d %e\n',left,En,size(S,1),Truncerr)
    V = conj(Vconj);
    % Renormalization
    LeftBlocks(end) = [];
    RightBlocks{end+1} = permute(tensorprod(tensorprod(V',newRightBlock,2,1),V,2,1),[1,3,2]);
  end

  % left end --> center
  while length(LeftBlocks) < L/2
    left = length(LeftBlocks);
    dimLeft = size(LeftBlocks{end},1);
    dimRight = size(RightBlocks{end},1);
    % form [[left]*]
    newLeftBlock = reshape(permute(tensorprod(LeftBlocks{end},HamMPO{left+1},3,3,"NumDimensionsA",3),...
                                   [1,3,2,4,5]), [dimLeft*d,dimLeft*d, size(HamMPO{L_act/2},4)]);
    % form [*[right]]
    newRightBlock = reshape(permute(tensorprod(HamMPO{left+2},RightBlocks{end},4,3,"NumDimensionsA",4),...
                                   [1,4,2,5,3]), [dimRight*d,dimRight*d, size(HamMPO{L_act/2+1},3)]);
    % Diagonalization
    dimSB = size(newLeftBlock,1)*size(newRightBlock,1);
    myHdotPsi = @(x) internal_HdotPsi(x,newLeftBlock,newRightBlock);
    [EigVec,En] = eigs(myHdotPsi,dimSB,[],1,'smallestreal');
    EigVec = reshape(EigVec,[dimLeft*d,dimRight*d]);  
    % Schmidt decomposition and truncation
    [U,S,Vconj] = svd(EigVec,'econ');
    Truncerr = 0;
    if size(S,1) > M
      U = U(:,1:M);
      Vconj = Vconj(:,1:M);
      S = diag(S);
      Truncerr = S(M+1:end)'*S(M+1:end);
      S = S(1:M);
    end
    fprintf('%d %d %d %e\n',left,En,size(S,1),Truncerr)
    V = conj(Vconj);
    % Renormalization
    LeftBlocks{end+1} = permute(tensorprod(tensorprod(U',newLeftBlock,2,1),U,2,1),[1,3,2]);
    RightBlocks(end) = [];
  end

end
plot(Energies)





function Y = internal_HdotPsi(X,newLeftBlock,newRightBlock)
  X = reshape(X,[size(newLeftBlock,1),size(newRightBlock,1)]);
  Y = tensorprod(newLeftBlock,X,2,1);
  Y = tensorprod(Y,newRightBlock,[2,3],[3,2]);
  Y = Y(:);   %reshape as a vector;
end
