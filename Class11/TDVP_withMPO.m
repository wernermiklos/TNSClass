L = 50;
J = 1;
Delta = 1;
dt = 0.01;
tmax = 1.0;
d = 2;   %local dimension
M = 32;
HamMPO = MPO_XXZ_full(L,J,Delta);
Nsweep = tmax / dt;

% Initialization  (here no chain-growth part)   
% we build a fully right canonical MPS. Init state: Neél state
MPS = struct();
MPS.LeftMatrices = {};
MPS.RightMatrices = {};
MPS.CoreMatrix = 0;
for i = 1:L
  V = zeros([1,2,1]);    %left_block, site_block, right_block
  if mod(i,2)    %up
    V(1,1,1) = 1;        %left_block, site_block, right_block
  else           %down
    V(1,2,1) = 1;
  end
  if i ~= L
    MPS.RightMatrices{end+1} = V;
  else
    MPS.CoreMatrix = V;
  end
end
% Initial renormalization (Block formation)
LeftBlocks = {1};    %dummy left block: 1x1x1 unit matrix
RightBlocks = {1};   %dummy right block: 1x1x1 unit matrix
for i = 1:(L-2)
  myMPOmatrix = HamMPO{end-i+1};
  tmp = tensorprod(RightBlocks{end},conj(MPS.RightMatrices{i}),1,3,'NumDimensionsA',3);
  tmp = tensorprod(tmp,myMPOmatrix,[4,2],[1,4],'NumDimensionsA',4);
  RightBlocks{end+1} = permute(tensorprod(tmp,MPS.RightMatrices{i},[1,3],[3,2],'NumDimensionsA',4),[1,3,2]);
end



Szdata = zeros(0,L);
times = zeros(0,1);
t = 0;
disp('SWEEPING')
disp('t left Sz_left Sz_left+1 M TrErr')
for sweep_iter = 1:Nsweep
  t = t + dt / 2;
  %center --> right end
  while length(MPS.RightMatrices) > 0

    % TWO SITE FORWARD
    left = length(LeftBlocks);
    dimLeft = size(LeftBlocks{end},1);
    dimRight = size(RightBlocks{end},1);
    % form [[left]*]
    newLeftBlock = reshape(permute(tensorprod(LeftBlocks{end},HamMPO{left},3,3,"NumDimensionsA",3),...
                                   [1,3,2,4,5]), [dimLeft*d,dimLeft*d, size(HamMPO{left},4)]);
    % form [*[right]]
    newRightBlock = reshape(permute(tensorprod(HamMPO{left+1},RightBlocks{end},4,3,"NumDimensionsA",4),...
                                   [1,4,2,5,3]), [dimRight*d,dimRight*d, size(HamMPO{left+1},3)]);
    % initvec
    Psi0 = reshape(tensorprod(MPS.CoreMatrix,MPS.RightMatrices{end},3,1),[dimLeft*d*dimRight*d,1]);

    
    % Diagonalization
    myHdotPsi = @(x) internal_HdotPsi(x,newLeftBlock,newRightBlock);
    Psi = expevolv_shifttrick(Psi0,+dt/2, myHdotPsi, NaN);   %NaN: adaptive exponential order

    Psi = reshape(Psi,[dimLeft*d,dimRight*d]);  
    % Schmidt decomposition and truncation
    [U,S,Vdag] = svd(Psi,'econ');
    Truncerr = 0;
    if size(S,1) > M
      U = U(:,1:M);
      Vdag = Vdag(:,1:M);
      S = diag(S);
      norm_orig = norm(S);
      Truncerr = S(M+1:end)'*S(M+1:end);
      S = S(1:M);
      S = S*norm_orig/norm(S);  %to keep norm
      S = diag(S);
    end
    %Measurement
    Psi = reshape(Psi,[dimLeft,d,d,dimRight]);
    rho1 = tensorprod(conj(Psi),Psi,[1,3,4],[1,3,4]);   %reduced DM of site 1
    rho2 = tensorprod(conj(Psi),Psi,[1,2,4],[1,2,4]);   %reduced DM of site 2
    Sz1 = trace(rho1*[1 0; 0 -1]);
    Sz2 = trace(rho2*[1 0; 0 -1]);
    fprintf('%d %d %d %d %d %e \n',t, left, Sz1, Sz2, size(S,1),Truncerr)
    % Renormalization
    MPS.LeftMatrices{end+1} = reshape(U,[dimLeft,d,size(U,2)]);
    MPS.CoreMatrix = reshape(S*Vdag',[size(S,1),d,dimRight]);
    MPS.RightMatrices(end) = [];
    if ~isempty(MPS.RightMatrices)
      LeftBlocks{end+1} = permute(tensorprod(tensorprod(U',newLeftBlock,2,1),U,2,1),[1,3,2]);
      % ONE SITE BACKWARD
      newLeftBlock = LeftBlocks{end};
      Psi0 = reshape(MPS.CoreMatrix,[size(MPS.CoreMatrix,1)*d*dimRight,1]);
      myHdotPsi = @(x) internal_HdotPsi(x,newLeftBlock,newRightBlock);
      Psi = expevolv_shifttrick(Psi0,-dt/2, myHdotPsi, NaN);   %NaN: adaptive exponential order
      MPS.CoreMatrix = reshape(Psi,[size(MPS.CoreMatrix,1),d,dimRight]);
      RightBlocks(end) = [];
    end
  end

  t = t + dt / 2;
  Szdata(end+1,:) = 0;
  times(end+1,1) = t;
  % Right to left sweep
  while length(MPS.LeftMatrices) > 0

    % TWO SITE FORWARD
    left = length(LeftBlocks);
    dimLeft = size(LeftBlocks{end},1);
    dimRight = size(RightBlocks{end},1);
    % form [[left]*]
    newLeftBlock = reshape(permute(tensorprod(LeftBlocks{end},HamMPO{left},3,3,"NumDimensionsA",3),...
                                   [1,3,2,4,5]), [dimLeft*d,dimLeft*d, size(HamMPO{left},4)]);
    % form [*[right]]
    newRightBlock = reshape(permute(tensorprod(HamMPO{left+1},RightBlocks{end},4,3,"NumDimensionsA",4),...
                                   [1,4,2,5,3]), [dimRight*d,dimRight*d, size(HamMPO{left+1},3)]);
    % initvec
    Psi0 = reshape(tensorprod(MPS.LeftMatrices{end},MPS.CoreMatrix,3,1),[dimLeft*d*dimRight*d,1]);

    
    % Diagonalization
    myHdotPsi = @(x) internal_HdotPsi(x,newLeftBlock,newRightBlock);
    Psi = expevolv_shifttrick(Psi0,+dt/2, myHdotPsi, NaN);   %NaN: adaptive exponential order

    Psi = reshape(Psi,[dimLeft*d,dimRight*d]);  
    % Schmidt decomposition and truncation
    [U,S,Vdag] = svd(Psi,'econ');
    Truncerr = 0;
    if size(S,1) > M
      U = U(:,1:M);
      Vdag = Vdag(:,1:M);
      S = diag(S);
      norm_orig = norm(S);
      Truncerr = S(M+1:end)'*S(M+1:end);
      S = S(1:M);
      S = S*norm_orig/norm(S);  %to keep norm
      S = diag(S);
    end
    %Measurement
    Psi = reshape(Psi,[dimLeft,d,d,dimRight]);
    rho1 = tensorprod(conj(Psi),Psi,[1,3,4],[1,3,4]);   %reduced DM of site 1
    rho2 = tensorprod(conj(Psi),Psi,[1,2,4],[1,2,4]);   %reduced DM of site 2
    Sz1 = trace(rho1*[1 0; 0 -1]);
    Sz2 = trace(rho2*[1 0; 0 -1]);
    Szdata(end,left) = Sz1;
    Szdata(end,left+1) = Sz2;
    fprintf('%d %d %d %d %d %e\n',t, left, Sz1, Sz2, size(S,1),Truncerr)
    % Renormalization
    MPS.RightMatrices{end+1} = reshape(Vdag',[size(Vdag,2),d,dimRight]);
    MPS.CoreMatrix = reshape(U*S,[dimLeft,d,size(S,2)]);
    MPS.LeftMatrices(end) = [];
    if ~isempty(MPS.LeftMatrices)
      RightBlocks{end+1} = permute(tensorprod(tensorprod(Vdag.',newRightBlock,2,1),conj(Vdag),2,1),[1,3,2]);
      % ONE SITE BACKWARD
      newRightBlock = RightBlocks{end};
      Psi0 = reshape(MPS.CoreMatrix,[dimLeft*d*size(MPS.CoreMatrix,3),1]);
      myHdotPsi = @(x) internal_HdotPsi(x,newLeftBlock,newRightBlock);
      Psi = expevolv_shifttrick(Psi0,-dt/2, myHdotPsi, NaN);   %NaN: adaptive exponential order
      MPS.CoreMatrix = reshape(Psi,[dimLeft,d,size(MPS.CoreMatrix,3)]);
      LeftBlocks(end) = [];
    end
  end

end
plot(times,Szdata(:,25));
hold on
plot(times,Szdata(:,26));





function Y = internal_HdotPsi(X,newLeftBlock,newRightBlock)
  X = reshape(X,[size(newLeftBlock,1),size(newRightBlock,1)]);
  Y = tensorprod(newLeftBlock,X,2,1);
  Y = tensorprod(Y,newRightBlock,[2,3],[3,2]);
  Y = Y(:);   %reshape as a vector;
end
