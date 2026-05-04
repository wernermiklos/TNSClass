%% iTEBD for the spin-1/2 XXZ chain in Vidal's canonical form
%  (version using tensorprod instead of bsxfun/reshape)
%
% Hamiltonian:  H = sum_i  Jxy * (Sx_i Sx_{i+1} + Sy_i Sy_{i+1})
%                         + Jz  *  Sz_i Sz_{i+1}
%
% Two-site unit cell (A,B), Vidal form:
%   ... la_B  Ga_A  la_A  Ga_B  la_B  Ga_A  la_A  Ga_B  la_B ...
%
% Index conventions (same as before):
%   Ga_A(s, a, b):   s = physical, a = left bond (dim chi_B), b = right bond (dim chi_A)
%   Ga_B(s, a, b):   s = physical, a = left bond (dim chi_A), b = right bond (dim chi_B)
%   la_A(a), la_B(a): Schmidt values (stored as vectors)
%
% First-order Trotter:  U(dt) = U_AB(dt) * U_BA(dt)
% Initial state: z-Neel  |up,down,up,down,...>  (product state, chi = 1)

clear; close all; clc;

%% -------------------- Parameters --------------------
Jxy    = 1.0;
Jz     = 1.0;
dt     = 0.05;
Nstep  = 400;
chiMax = 256;
tol    = 1e-12;

d = 2;

%% -------------------- Local operators --------------------
Sx = 0.5*[0 1; 1 0];
Sy = 0.5*[0 -1i; 1i 0];
Sz = 0.5*[1 0; 0 -1];

% Two-site Hamiltonian on a bond, basis |s1 s2>
h2 = Jxy*( kron(Sx,Sx) + kron(Sy,Sy) ) + Jz*kron(Sz,Sz);

% Trotter gate U(s1', s2', s1, s2) = <s1' s2'| exp(-i h dt) |s1 s2>
U = reshape( expm(-1i*dt*h2), [d d d d] );

%% -------------------- Initial state: z-Neel --------------------
GA = zeros(d,1,1);  GA(1,1,1) = 1;     % |up>
GB = zeros(d,1,1);  GB(2,1,1) = 1;     % |down>
lA = 1;
lB = 1;

%% -------------------- Storage --------------------
tlist   = (1:Nstep)*dt;
SzA_arr = zeros(1,Nstep);
SzB_arr = zeros(1,Nstep);
SvnA_arr = zeros(1,Nstep);   % entropy across the A-bond
SvnB_arr = zeros(1,Nstep);   % entropy across the B-bond

%% -------------------- Time evolution --------------------
for n = 1:Nstep
  disp(n);
    % Gate on AB bond (outer env lB, inner bond lA)
    [GA, lA, GB] = applyGate(GA, lA, GB, lB, U, chiMax, tol);
    % Gate on BA bond (outer env lA, inner bond lB)
    [GB, lB, GA] = applyGate(GB, lB, GA, lA, U, chiMax, tol);

    % Measure <S^z> on each sublattice
    SzA_arr(n) = measureOneSite(GA, lB, lA, Sz);
    SzB_arr(n) = measureOneSite(GB, lA, lB, Sz);

    % Von Neumann entropy from Schmidt values:  S = -sum p_i log(p_i), p_i = lambda_i^2
    SvnA_arr(n) = vonNeumannEntropy(lA);
    SvnB_arr(n) = vonNeumannEntropy(lB);
end

%% -------------------- Plot --------------------
figure;
plot(tlist, SzA_arr, 'b-', 'LineWidth', 1.4); hold on;
plot(tlist, SzB_arr, 'r-', 'LineWidth', 1.4);
xlabel('time  t');
ylabel('\langle S^z \rangle');
legend('sublattice A', 'sublattice B', 'Location','best');
title(sprintf('iTEBD, XXZ chain, J_{xy}=%.2f, J_z=%.2f, \\chi_{max}=%d, dt=%.3f', ...
               Jxy, Jz, chiMax, dt));
grid on;

figure;
plot(tlist, SvnA_arr, 'b-', 'LineWidth', 1.4); hold on;
plot(tlist, SvnB_arr, 'r--', 'LineWidth', 1.4);
xlabel('time  t');
ylabel('von Neumann entropy  S');
legend('A-bond', 'B-bond', 'Location','best');
title(sprintf('Entanglement entropy across a bond (\\chi_{max}=%d)', chiMax));
grid on;

save(['Data_chi',num2str(chiMax),'.mat'],'tlist',"SvnA_arr")

%% ================================================================
%% ===================  Helper functions  =========================
%% ================================================================

function [G1new, lnew, G2new] = applyGate(G1, l1, G2, l2, U, chiMax, tol)
% Apply a two-site gate U on the bond between sites 1 and 2 in Vidal form.
%
% Environment:  ... l2 -- G1 -- l1 -- G2 -- l2 ...
%   G1(s1, a, b):  a of dim chL = length(l2),  b of dim chM = length(l1)
%   G2(s2, b, c):  b of dim chM,               c of dim chR = length(l2)
%   U(s1', s2', s1, s2)

    d   = size(G1,1);
    chL = size(G1,2);
    chR = size(G2,3);

    % ---- 1. Absorb Schmidt values:  T1 = l2 * G1 * l1,   T2 = G2 * l2 ----
    % Multiply a diagonal matrix on a bond == contract with diag(l).
    % Use tensorprod with a diagonal: tensorprod(diag(l), X, 2, axis).
    T1 = tensorprod( diag(l2), G1, 2, 2 );   % contracts l2's col-index with G1's left bond
    % T1 now has indices (a, s1, b).  Bring it to (s1, a, b):
    T1 = permute(T1, [2 1 3]);
    T1 = tensorprod( T1, diag(l1), 3, 1 );   % contract right bond of T1 with l1's row-index
    % T1 indices: (s1, a, b)

    T2 = tensorprod( G2, diag(l2), 3, 1 );   % (s2, b, c)

    % ---- 2. Build theta by contracting on the middle bond b ----
    % theta(s1, a, s2, c) = sum_b T1(s1,a,b) * T2(s2,b,c)
    theta = tensorprod(T1, T2, 3, 2);        % indices (s1, a, s2, c)

    % ---- 3. Apply gate: contract (s1,s2) of theta with (s1,s2) of U ----
    % U has indices (s1', s2', s1, s2); contract U's last two with theta's (s1,s2).
    % After contraction, leftover indices of U are (s1', s2'); of theta are (a, c).
    % tensorprod result ordering: [leftover(U), leftover(theta)] = (s1', s2', a, c)
    theta = tensorprod(U, theta, [3 4], [1 3]);   % (s1', s2', a, c)

    % ---- 4. SVD: group (s1', a) vs (s2', c) ----
    % Reshape is unavoidable here because svd works on matrices.
    theta = permute(theta, [1 3 2 4]);       % (s1', a, s2', c)
    M     = reshape(theta, [d*chL, d*chR]);

    [Usvd, S, Vsvd] = svd(M, 'econ');
    svals = diag(S);

    chi_new = max(1, min(chiMax, sum(svals > tol)));
    Usvd  = Usvd(:, 1:chi_new);
    Vsvd  = Vsvd(:, 1:chi_new);
    svals = svals(1:chi_new);

    svals = svals / norm(svals);
    lnew  = svals(:).';

    % ---- 5. Rebuild Gammas ----
    % Usvd is (s1'*a, chi_new) -> A(s1', a, chi_new)
    A = reshape(Usvd, [d, chL, chi_new]);
    % Remove outer l2 on the left:  G1new = l2^{-1} * A
    invL = 1 ./ l2;  invL(~isfinite(invL)) = 0;
    G1new = tensorprod( diag(invL), A, 2, 2 );   % (a, s1', chi_new)
    G1new = permute(G1new, [2 1 3]);             % (s1', a, chi_new)

    % Vsvd is (s2'*c, chi_new) -> V^dagger has shape (chi_new, s2', c)
    B = reshape(Vsvd', [chi_new, d, chR]);
    B = permute(B, [2 1 3]);                     % (s2', chi_new, c)
    % Remove outer l2 on the right
    invR = 1 ./ l2;  invR(~isfinite(invR)) = 0;
    G2new = tensorprod( B, diag(invR), 3, 1 );   % (s2', chi_new, c)
end


function val = measureOneSite(G, lL, lR, Op)
% <Op> at a site with Vidal tensor G(s,a,b) and Schmidt environments
% lL on the left bond (a) and lR on the right bond (b).
%
% rho(s,s') = sum_{a,b} lL(a)^2 * G(s,a,b) * conj(G(s',a,b)) * lR(b)^2

    % Weight G by lL (left) and lR (right):
    GW = tensorprod( diag(lL), G,  2, 2 );   % (a, s, b)
    GW = permute(GW, [2 1 3]);               % (s, a, b)
    GW = tensorprod( GW, diag(lR), 3, 1 );   % (s, a, b)

    % rho(s, s') = sum_{a,b} GW(s,a,b) * conj(GW(s',a,b))
    rho = tensorprod( GW, conj(GW), [2 3], [2 3] );   % (s, s')

    val = real( trace(Op * rho) );
end


function S = vonNeumannEntropy(lambda)
% Entanglement entropy from Schmidt values:
%   S = -sum_i p_i log(p_i),   p_i = lambda_i^2
% Assumes lambda is normalized (sum of squares = 1).

    p = lambda(:).^2;
    p = p(p > 0);            % drop zeros to avoid 0*log(0)
    S = -sum( p .* log(p) );
end