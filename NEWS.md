# KPC 0.1.3

## Changes in this version:
- Added the GitHub URL to the `DESCRIPTION` file.
- Enabled parallel support for `KFOCI()` and `KPCRKHS_VS()` on Windows.
- Added functionality to pre-condition on a set of variables in `KFOCI()`.

---

# KPC 0.1.2

## Changes in this version:
- Added the URL to the published paper in the `DESCRIPTION` file.
- Updated the default choice of `K` for the K-nearest neighbor (K-NN) graph in the variable selection algorithm `KFOCI()` to `0.05n`. Additional details and suggestions for choosing `K` have been included in the documentation.
- Enabled parallel execution by default for `KFOCI()` and `KPCRKHS_VS()`.
