function [K, P, nH, m, ns, nz, rhs] = assembleK3p(itnumber, varargin)
  % Assemble
  %       [ H + rho*I   -I      J'   ] } n
  % K3p = [     Z        X           ] } ns
  %       [     J           -delta I ] } m
  %
  % and return the corresponding right-hand side and a preconditioner P.
  %
  % `[K, P, rhs] = assembleK3p(itnumber)`
  %   returns P = I.
  %
  % `[K, P, rhs] = assembleK3p(itnumber, precon, arg1, ..., argN)`
  %   calls `precon(rho, delta, H, J, Z, X, arg1, ..., argN). The
  %   function `precon()` should return an appropriate preconditioner
  %   for K3p.
  %
  [rho, delta, H, J, Z, X, rhs] = read_blocks(itnumber);

  % I is really a rectangular block shaped like Z'.
  I = Z';
  I(find(I > 0)) = 1.0;

  [m, n] = size(J);
  ns = size(Z, 1);      % Number of slack variables.
  nH = size(H, 1);
  n = nH - ns;  % Number of original variables.

  K = [ H + tril(H,-1)'    -I               J'               ; ...
        Z.^2                X               sparse(ns, m)    ; ...
        J                   sparse(m, ns)  -delta * speye(m)];

  % Unscale right-hand side.
  rhs(end-ns+1:end) = Z(:, n+1:end) * rhs(end-ns+1:end);

  % Reorder right-hand side.
  rhs1 = rhs(1:n+ns);         % length n+ns
  rhs2 = rhs(n+ns+1:n+ns+m);  % length m
  rhs3 = rhs(n+ns+m+1:end);   % length ns
  rhs = [rhs1 ; rhs3 ; rhs2];

  % Obtain preconditioner.
  if nargin > 1
    [P, nz] = varargin{1}(rho, delta, H, J, Z, X, varargin{2:end});
  else
    P = speye(size(K));
    nz = size(K, 1);
  end
end
