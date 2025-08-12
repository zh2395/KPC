# KPC 0.1.3

## Changes in this version:
- Added the GitHub URL to the `DESCRIPTION` file.
- Enabled parallel support for `KFOCI()` and `KPCRKHS_VS()` on Windows.
- Added functionality to pre-condition on a set of variables in `KFOCI()`.

---

# KPC 0.1.2

## Changes in this version:
- The URL to the published paper is added to the `DESCRIPTION`.
- The default choice of K for the K-nearest neighbor (K-NN) graph in the variable selection algorithm `KFOCI()` is changed to 0.05n. More details of our suggestions on the choice of K are added to the documentation.
- `KFOCI()` and `KPCRKHS_VS()` will run in parallel by default.
