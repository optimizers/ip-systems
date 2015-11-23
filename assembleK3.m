function [K, P, nz, rhs] = assembleK3(itnumber, varargin)
  % Assemble
  %      [ H + rho*I     J'    -I ]
  % K3 = [     J     -delta*I     ]
  %      [     Z                X ]
  %
  % and return the corresponding right-hand side and a preconditioner P.
  %
  % `[K, P, rhs] = assembleK3(itnumber)`
  %   returns P = I.
  %
  % `[K, P, rhs] = assembleK3(itnumber, precon, arg1, ..., argN)`
  %   calls `precon(rho, delta, H, J, Z, X, arg1, ..., argN). The
  %   function `precon()` should return an appropriate preconditioner
  %   for K3.
  %
  [rho, delta, H, J, Z, X, rhs] = read_blocks(itnumber);

  % I is really a rectangular block shaped like Z'.
  I = Z';
  I(find(I > 0)) = 1.0;

  [m, n] = size(J);
  ns = size(Z, 1);      % Number of slack variables.
  n = size(H, 1) - ns;  % Number of original variables.

  K = [ H + tril(H,-1)'      J'                -I              ; ...
        J                   -delta * speye(m)   sparse(m, ns)  ; ...
        Z.^2                 sparse(ns, m)      X             ];

  % Unscale right-hand side.
  rhs(end-ns+1:end) = Z(:, n+1:end) * rhs(end-ns+1:end);

  % Obtain preconditioner.
  if nargin > 1
    [P, nz] = varargin{1}(rho, delta, H, J, Z, X, varargin{2:end});
  else
    P = speye(size(K));
    nz = size(K, 1);
  end
end
