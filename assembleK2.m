function [K, P, nz, rhs] = assembleK2(itnumber, varargin)
  % Assemble
  %      [ H + rho*I + X^{-1} Z     J'     ]
  % K2 = [         J              -delta*I ]
  %
  % and return the corresponding right-hand side and a preconditioner P.
  %
  % `[K, P, rhs] = assembleK2(itnumber)`
  %   returns P = I.
  %
  % `[K, P, rhs] = assembleK2(itnumber, precon, arg1, ..., argN)`
  %   calls `precon(rho, delta, H, J, Z, X, arg1, ..., argN). The
  %   function `precon()` should return an appropriate preconditioner
  %   for K2.
  %
  [rho, delta, H, J, Z, X, rhs] = read_blocks(itnumber);

  [m, n] = size(J);
  ns = size(Z, 1);      % Number of slack variables.
  n = size(H, 1) - ns;  % Number of original variables.

  XinvZ = Z' * (X \ Z) + rho * speye(size(Z, 2));
  K = [ H + tril(H,-1)' + XinvZ      J'                ; ...
        J                           -delta * speye(m) ];

  % Reduce right-hand side.
  rhs(1:n+ns) = rhs(1:n+ns) - Z' * (X \ rhs(n+ns+m+1:n+ns+m+ns));
  rhs = rhs(1:n+ns+m);

  % Obtain preconditioner.
  if nargin > 1
    [P, nz] = varargin{1}(rho, delta, H, J, Z, X, varargin{2:end});
  else
    P = speye(size(K));
    nz = size(K, 1);
  end
end
