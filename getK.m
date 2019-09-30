function [K, P, n, m, ns, nz, rhs] = getK(assembler, problem, itnumber, varargin)
  % Usage:
  %   [K, P, n, m, ns, nz, rhs] = getK(@assembleK3, 'cvxqp1_s', 5)
  %   [K, P, n, m, ns, nz, rhs] = getK(@assembleK3, 'cvxqp1_s', 5, @my_preconditioner, args...)
  %
  cur = pwd;
  cd(fullfile(cur, 'data', sprintf('%s', problem), '3x3', ...
              sprintf('iter_%d', itnumber)));

  [K, P, n, m, ns, nz, rhs] = assembler(itnumber, varargin{:});
  cd(cur);
end
