# Interior-Point Systems

This repository contains a collection of linear systems arising from
interior-point methods for quadratic optimization in
[MatrixMarket](http://math.nist.gov/MatrixMarket/formats.html#MMformat) format.
The distinguishing features of the collection are that

* systems contain accompanying right-hand sides
* systems are supplied in the form of their blocks, allowing users to solve
  several equivalent formulations of the same systems
* sets of related systems are supplied, generated during the iterations of an
  interior-point method applied to the same optimization problem

## Usage

A Matlab interface to the systems is provided. From the top-level folder,
```matlab
[P, K, nz, rhs] = getK(@my_assembler, problem, iter, @my_preconditioner, args...)
```

returns the system generated from problem `problem` at interior-point iteration
`iter` in `K` and `rhs`. The assembler `@my_assembler` assembles the system
from its blocks, as read by `read_blocks()`. Example assemblers are provided in
[`assembleK3`](https://github.com/optimizers/ip-systems/blob/master/assembleK3.m),
[`assembleK35`](https://github.com/optimizers/ip-systems/blob/master/assembleK35.m)
and
[`assembleK2`](https://github.com/optimizers/ip-systems/blob/master/assembleK2.m).
If supplied, `my_preconditioner` should return an adequate preconditioner `P`
together with a measure of its "complexity" (e.g., its number of nonzeros) in
`nz`. Additional arguments `args...` are passed unchanged to
`my_preconditioner`. If no preconditioner is supplied, `P` is set to a sparse
identity matrix.

See the technical report below for examples.

## Citing this Collection

If you use this collection in your research, please cite the following sources

1. Orban D., A Collection of Linear Systems Arising from Interior-Point Methods
   for Quadratic Optimization, Cahier du GERAD G-2015-00, GERAD, Montreal,
   Canada, 2015.
   [BibTeX](https://github.com/optimizers/ip-systems/blob/master/misc/citing.bib#L1).
2. Orban D., A Collection of Linear Systems Arising from Interior-Point Methods
   for Quadratic Optimization, online data set, 2015.
   [BibTeX](https://github.com/optimizers/ip-systems/blob/master/misc/citing.bib#L12).
