% Example script to cycle through all problems.
% clear all;

problems = { ...
  'aug2d', ...
  'aug2dc', ...
  'aug2dcqp', ...
  'aug2dqp', ...
  'aug3d', ...
  'aug3dc', ...
  'aug3dcqp', ...
  'aug3dqp', ...
  'cvxqp1_s', ...
  'cvxqp1_m', ...
  'cvxqp1_l', ...
  'cvxqp2_s', ...
  'cvxqp2_m', ...
  'cvxqp2_l', ...
  'cvxqp3_s', ...
  'cvxqp3_m', ...
  'cvxqp3_l', ...
  'dual1', ...
  'dual2', ...
  'dual3', ...
  'dual4', ...
  'dualc1', ...
  'dualc2', ...
  'dualc5', ...
  'dualc8', ...
  'genhs28', ...
  'gouldqp2', ...
  'gouldqp3', ...
  'hs118', ...
  'hs21', ...
  'hs21mod', ...
  'hs268', ...
  'hs35', ...
  'hs35mod', ...
  'hs51', ...
  'hs52', ...
  'hs53', ...
  'hs76', ...
  'hues-mod', ...
  'huestis', ...
  'ksip', ...
  'liswet1', ...
  'liswet10', ...
  'liswet11', ...
  'liswet12', ...
  'liswet2', ...
  'liswet3', ...
  'liswet4', ...
  'liswet5', ...
  'liswet6', ...
  'liswet7', ...
  'liswet8', ...
  'liswet9', ...
  'lotschd', ...
  'mosarqp1', ...
  'mosarqp2', ...
  'powell20', ...
  'primal1', ...
  'primal2', ...
  'primal3', ...
  'primal4', ...
  'primalc1', ...
  'primalc2', ...
  'primalc5', ...
  'primalc8', ...
  'qpcblend', ...
  'qpcboei1', ...
  'qpcboei2', ...
  'qpcstair', ...
  's268', ...
  'stcqp1', ...
  'stcqp2', ...
  'tame', ...
  'ubh1', ...
  'yao', ...
  'zecevic2', ...
};

assemblers = {  ...
  @assembleK3,  ...
  @assembleK35, ...
  @assembleK2,  ...
};

fprintf('%14s  %23s  %23s  %23s\n', '', 'K3', 'K3.5', 'K2');
fprintf('%10s  %2s  %5s  %7s  %7s  %5s  %7s  %7s  %5s  %7s  %7s\n', ...
        'problem', 'it', ...
        'size', 'dens(%)', 'cond', ...
        'size', 'dens(%)', 'cond', ...
        'size', 'dens(%)', 'cond');

ip_iter = 10;  % interior-point iteration.

for p = 1 : length(problems)
  problem = problems{p};
  
  iter_subdir = fullfile('data', sprintf('%s', problem), '3x3', ...
                         sprintf('iter_%d', ip_iter));
  
  if exist(iter_subdir, 'dir') == 7
    fprintf('%10s  %2d', problem, ip_iter);

    for a = 1 : length(assemblers)
      assembler = assemblers{a};
      [K, P, nz, rhs] = getK(assembler, problem, ip_iter);

      n = size(K, 1);
      dens = 100 * nnz(K) / (n * (n+1) / 2);
      fprintf('  %5d  %7.1e  %7.1e', n, dens, condest(K));
    end
    fprintf('\n');
  end
end