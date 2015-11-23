function [K, P, nz, rhs] = assembleK35(itnumber, varargin)
  % Assemble
  %        [ H + rho*I     J'    -Z^{1/2} ]
  % K3.5 = [     J     -delta*I           ]
  %        [ -Z^{1/2}               -X    ]
  %
  % and return the corresponding right-hand side and a preconditioner P.
  %
  % `[K, P, rhs] = assembleK35(itnumber)`
  %   returns P = I.
  %
  % `[K, P, rhs] = assembleK35(itnumber, precon, arg1, ..., argN)`
  %   calls `precon(rho, delta, H, J, Z, X, arg1, ..., argN). The
  %   function `precon()` should return an appropriate preconditioner
  %   for K3.5.
  %
  [rho, delta, H, J, Z, X, rhs] = read_blocks(itnumber);

  [m, n] = size(J);
  ns = size(Z, 1);      % Number of slack variables.

  K = [ H + tril(H,-1)'      J'                -Z'             ; ...
        J                   -delta * speye(m)   sparse(m, ns)  ; ...
       -Z                    sparse(ns, m)     -X             ];

  % Obtain preconditioner.
  if nargin > 1
    [P, nz] = varargin{1}(rho, delta, H, J, Z, X, varargin{2:end});
  else
    P = speye(size(K));
    nz = size(K, 1);
  end
end
