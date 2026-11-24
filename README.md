# DisorderBHZ

## Reduced spectral post-processing

The `spectrum_Abhz_2d` driver rebuilds one seed's free or converged
mean-field Hamiltonian and computes local diagonal spectra with LAPACK
or Chebychev. It never writes the full
`Nlat x Lfreq x 4` site-frequency array. For each frequency and local
channel it writes only Welford sufficient statistics and numerical
diagnostics directly under

```text
PREFIX/seed_<IDUM>/
```

The explicit channel names are `up_orb1`, `up_orb2`, `down_orb1`,
`down_orb2`, and `total`. The `total` spectrum is formed locally before
its logarithm, so its geometric average is well defined.

The spectral controls are independent of the Chebychev Fermi-projector
controls:

```text
SPECTRUM_SOURCE=free              ! free or mean_field
SPECTRUM_METHOD=lapack            ! lapack or chebychev
SPECTRUM_ORDER=1024               ! Chebychev resolvent order
SPECTRUM_RHO_FLOOR=1d-300         ! used only inside log(rho)
SPECTRUM_NEGATIVE_TOLERANCE=1d-12
```

`SPECTRUM_SOURCE=mean_field` reads `params.restart` literally and
requires the last error in `mf_history.dat` to be smaller than
`IT_ERROR`. It does not add `SB_FIELD` again. The disorder vector saved
as `erandom_<IDUM>.used` is reused as the authoritative realization,
with `.restart` retaining precedence when both exist.

Build and run, using the same input file and command-line override
conventions as the solution drivers:

```text
make spectrum
spectrum_Abhz_2d INPUTFILE=inputABHZ.conf \
  SPECTRUM_SOURCE=free SPECTRUM_METHOD=lapack
```

The `spectrum` target also installs the executable Python merger as
`$(DIREXE)/merge_spectrum_stats.py`; `DIREXE` defaults to `~/.bin`.

Merge compatible seed directories with equal seed weights:

```text
merge_spectrum_stats.py \
  run/seed_* -o run/spectrum_ensemble
```

The merger writes separate arithmetic, geometric, and diagnostic files
for every channel. It also creates one spectrum plot per channel and a
combined `rho_typical/rho_average` plot; use `--aggregate-only` to skip
figures. Spatial population standard deviations and seed-to-seed sample
standard deviations/SEM remain distinct. All table headers number their
columns explicitly.
