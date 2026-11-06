# Reference and experiment scope

S. A. Vorobyov, A. B. Gershman, and Z.-Q. Luo, “Robust Adaptive Beamforming Using Worst-Case Performance Optimization: A Solution to the Signal Mismatch Problem,” *IEEE Transactions on Signal Processing*, 51(2), 313–324, 2003. [DOI: 10.1109/TSP.2002.806865](https://doi.org/10.1109/TSP.2002.806865).

The signal look-direction mismatch is Example 2 in Section VI-B. The paper's example uses 3° presumed and 5° actual directions. The repository's main plot uses 2° presumed and 5° actual; the separate sweep uses 3° and 5°. The implementation is an adaptation of the paper's method. Its solver, uncertainty radius, baseline set and sweep ranges are described in the scripts.
