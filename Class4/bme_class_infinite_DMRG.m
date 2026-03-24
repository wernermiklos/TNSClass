%===========================================================
% Infinite-system DMRG for the 1D Transverse Field Ising Model
%===========================================================

% Parameters
J = 1.0;          % Ising coupling
h = 1.0;          % transverse field
M = 4;            % number of kept states
nSteps = 19;       % number of infinite-system growth steps

% Local operators
sx = sparse([0 1; 1 0]);
sz = sparse([1 0; 0 -1]);
id2 = speye(2);

%-----------------------------------------------------------
% Initial one-site block
% H_block = - h * sigma_x
% edge operator = sigma_z on the last site of the block
%-----------------------------------------------------------
H_block = -h * sx;
sz_edge = sz;

fprintf('Infinite chain DMRG for the transverse-field ising model\n');
fprintf('----------------------------------------\n');

for step = 1:nSteps

    %-------------------------------------------------------
    % 1. Enlarge system block by one site
    %
    % H_new = H_block \otimes I
    %       - J * sz_edge \otimes sz
    %       - h * I \otimes sx
    %-------------------------------------------------------
    dimBlock = size(H_block,1);

    H_sys = kron(H_block, id2) ...
          - J * kron(sz_edge, sz) ...
          - h * kron(eye(dimBlock), sx);

    % New edge operator of the enlarged system block
    sz_sys_edge = kron(eye(dimBlock), sz);

    %-------------------------------------------------------
    % 2. Symmetric infinite algorithm:
    %    environment block = mirror copy of system block
    %-------------------------------------------------------
    H_env = H_sys;
    sz_env_edge = sz_sys_edge;

    dimSys = size(H_sys,1);
    dimEnv = size(H_env,1);

    %-------------------------------------------------------
    % 3. Build superblock Hamiltonian
    %
    % H_SB = H_S \otimes I_E + I_S \otimes H_E
    %        - J * sz_edge(S) \otimes sz_edge(E)
    %-------------------------------------------------------
    H_super = kron(H_sys, eye(dimEnv)) ...
            + kron(eye(dimSys), H_env) ...
            - J * kron(sz_sys_edge, sz_env_edge);

    %-------------------------------------------------------
    % 4. Diagonalize superblock Hamiltonian
    %    Ground state = eigenvector with lowest eigenvalue
    %-------------------------------------------------------
    [Vsuper, Dsuper] = eigs(H_super);
    energies = diag(Dsuper);
    [E0, idx0] = min(energies);
    psi = Vsuper(:, idx0);

    %-------------------------------------------------------
    % 5. Reshape ground state into matrix Psi_{i,j}
    %    where i = system basis index, j = environment basis index
    %-------------------------------------------------------
    Psi = reshape(psi, [dimSys, dimEnv]);

    %-------------------------------------------------------
    % 6. Reduced density matrix of the system
    %
    % rho_S = Psi * Psi^\dagger
    %-------------------------------------------------------
    rho = Psi * Psi';

    % Ensure hermiticity numerically
    rho = 0.5 * (rho + rho');

    %-------------------------------------------------------
    % 7. Diagonalize rho and keep the M largest eigenvectors
    %
    % rho = U D U^\dagger
    %-------------------------------------------------------
    [U_rho, D_rho] = eig(rho);
    w = diag(D_rho);

    % Sort eigenvalues descending
    [w_sorted, order] = sort(w, 'descend');

    mKeep = min(M, length(w_sorted));
    keep = order(1:mKeep);

    % Truncation matrix U
    U = U_rho(:, keep);


    %-------------------------------------------------------
    % 8. Renormalize operators
    %
    % H_block <- U^\dagger H_sys U
    % sz_edge <- U^\dagger sz_sys_edge U
    %-------------------------------------------------------
    H_block = U' * H_sys * U;
    sz_edge = U' * sz_sys_edge * U;



    %-------------------------------------------------------
    % 9. Diagnostics
    %-------------------------------------------------------
    truncError = 1 - sum(w_sorted(1:mKeep));
    chainLength = 2 * (step + 1);   % starts from 2 sites, grows by 2

    fprintf('Step %2d | L = %2d | E0 = % .10f | discarded weight = %.3e\n', ...
            step, chainLength, E0, truncError);
end
