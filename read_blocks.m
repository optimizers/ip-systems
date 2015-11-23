function [rho, delta, H, J, Z, X, rhs] = read_blocks(itnumber)
  %
  % [rho, delta, H, J, Z, X, rhs] = read_blocks(itnumber)
  % Read blocks of a system at interior-point iteration `itnumber`.
  %
  % rho:   current value of primal regularization parameter
  % delta: current value of dual regularization parameter
  % H:     lower triangle of H + rho * I
  % J:     Jacobian of problem in slack form [J1 J2]
  % Z:     dual variables of problem in slack form [0 Z^{1/2}]
  % X:     diagonal matrix of slack variables
  % rhs:   right-hand side.

  rd = load(sprintf('rho_delta_%d.dat', itnumber));
  rho = rd(1); delta = rd(2);

  % This is the lower triangle of H only.
  % H already contains rho*I.
  H = mmread(sprintf('H+rhoI_%d.mtx', itnumber));
  H = -H;

  % Constraint matrix, including slack contributions.
  J = mmread(sprintf('J1J2_%d.mtx', itnumber));

  % This is sqrt(Z).
  Z = mmread(sprintf('Zsqrt_%d.mtx', itnumber));

  % Slack variables.
  X = mmread(sprintf('S_%d.mtx', itnumber));

  % This rhs really corresponds to K3.5.
  rhs = load(sprintf('rhs_%d.rhs', itnumber));
end
