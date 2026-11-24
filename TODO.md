# AndersonBHZ implementation roadmap

The current code provides the dense LAPACK reference path and the three
spin representations:

- `spin_blocks=1`: solve the up block and reconstruct the down block;
- `spin_blocks=2`: solve the up and down blocks independently;
- `spin_blocks=0`: solve the complete spin-orbital Hamiltonian.

The next development phase is the large-system Chebychev path.

## Chebychev projector and LCM

1. Define a matrix-free Hamiltonian action in `BHZ_SOLVE`. The driver
   must continue to own the physical Hamiltonian and pass either a sparse
   matrix or an equivalent matrix-vector callback to the solver.

2. Separate spectral-bound estimation from the projector recurrence.
   Retain the current sparse Gershgorin estimate as a safe reference and
   add a tighter optional estimate only if it gives a measured benefit.

3. Generalize the Chebychev Fermi projector so that `BHZ_LCM` can apply
   `P`, `PXP`, `PYP`, and the periodic position operators without
   constructing the dense occupied-state matrix `Uocc`.

4. Implement the three `spin_blocks` paths without adding separate public
   LCM interfaces. `spin_blocks=0` must treat the projected-spin operator
   `P Sz P` consistently and must check that its spectral gap is open.

5. Decide whether the `P Sz P` sector projectors should be obtained by a
   second polynomial filter or by a sparse eigensolver. This decision must
   be based on memory cost, number of matrix-vector products, and numerical
   stability near a closing projected-spin gap.

## Validation

1. Add small deterministic OBC and PBC regression cases for
   `spin_blocks=1,2,0`.

2. Compare Chebychev densities, Chern numbers, and local markers against
   the LAPACK reference while increasing the polynomial order.

3. Verify the four-commutator PBC expression independently and test the
   time-reversal relation between the two spin sectors.

4. Record peak memory and matrix-vector counts as functions of system size
   and polynomial order. Compilation alone is not considered validation.

5. Keep the mean-field loop explicit in `mf_Abhz_2d.f90`; only the linear
   solution step and projector operations belong in `BHZ_SOLVE`.
