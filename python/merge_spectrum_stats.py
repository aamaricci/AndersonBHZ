#!/usr/bin/env python3
"""Merge reduced Anderson-BHZ spectra from independent seed folders."""

import argparse
import math
import sys
from pathlib import Path

import numpy as np


CHANNELS = ("up_orb1", "up_orb2", "down_orb1", "down_orb2", "total")
CHANNEL_LABELS = {
    "up_orb1": r"$\uparrow$, orbital 1",
    "up_orb2": r"$\uparrow$, orbital 2",
    "down_orb1": r"$\downarrow$, orbital 1",
    "down_orb2": r"$\downarrow$, orbital 2",
    "total": "total",
}
CONSISTENT_KEYS = (
    "format", "source", "method", "nlat", "lfreq", "wmin", "wmax",
    "mu", "ieta", "spectrum_order", "rho_floor",
    "negative_tolerance", "spin_blocks", "channels",
)


def read_metadata(folder):
    """Read the whitespace-separated seed metadata contract."""
    path = folder / "spectrum_metadata.dat"
    values = {}
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: missing value")
            values[fields[0]] = " ".join(fields[1:])
    missing = set(CONSISTENT_KEYS + ("seed",)) - set(values)
    if missing:
        raise ValueError(f"{path}: missing metadata {sorted(missing)}")
    return values


def numeric(value):
    """Parse Fortran or Python floating-point text."""
    return float(value.replace("D", "E").replace("d", "e"))


def metadata_equal(left, right, key):
    """Compare numeric controls tightly and categorical controls exactly."""
    if key in {"wmin", "wmax", "mu", "ieta", "rho_floor",
               "negative_tolerance"}:
        return math.isclose(numeric(left), numeric(right),
                            rel_tol=1e-13, abs_tol=1e-15)
    return left == right


def check_metadata(folders):
    """Reject duplicate seeds and incompatible physical/numerical grids."""
    metadata = [read_metadata(folder) for folder in folders]
    reference = metadata[0]
    for folder, current in zip(folders[1:], metadata[1:]):
        for key in CONSISTENT_KEYS:
            if not metadata_equal(reference[key], current[key], key):
                raise ValueError(
                    f"{folder}: metadata {key}={current[key]!r} differs "
                    f"from {reference[key]!r}")
    seeds = [int(item["seed"]) for item in metadata]
    if len(seeds) != len(set(seeds)):
        raise ValueError(f"duplicate seed IDs: {seeds}")
    return reference, seeds


def read_channel(folder, channel, nlat, lfreq):
    """Read and validate one channel's sufficient statistics."""
    stats_path = folder / f"spectrum_{channel}_stats.dat"
    diag_path = folder / f"spectrum_{channel}_diagnostics.dat"
    stats = np.loadtxt(stats_path, ndmin=2)
    diag = np.loadtxt(diag_path, ndmin=2)
    if stats.shape != (lfreq, 6):
        raise ValueError(f"{stats_path}: expected {lfreq} rows and 6 columns")
    if diag.shape != (lfreq, 5):
        raise ValueError(f"{diag_path}: expected {lfreq} rows and 5 columns")
    if not np.isfinite(stats).all() or not np.isfinite(diag).all():
        raise ValueError(f"{folder}: nonfinite {channel} statistics")
    if (np.any(np.diff(stats[:, 0]) <= 0)
            or not np.array_equal(stats[:, 1], np.full(lfreq, nlat))
            or np.any(stats[:, 3] < 0)
            or np.any(stats[:, 5] < 0)
            or not np.allclose(stats[:, 0], diag[:, 0], rtol=0, atol=1e-13)):
        raise ValueError(f"{folder}: invalid {channel} statistics")
    integer_diag = diag[:, (1, 3)]
    if (np.any(integer_diag != np.rint(integer_diag))
            or np.any(integer_diag < 0)
            or np.any(integer_diag > nlat)):
        raise ValueError(f"{folder}: invalid {channel} diagnostic counts")
    return stats, diag


def merge_population(parts, mean_column, m2_column, count_column=1):
    """Merge vector Welford states with the stable parallel formula."""
    total = 0
    mean = np.zeros(parts[0].shape[0])
    m2 = np.zeros_like(mean)
    for part in parts:
        count = int(part[0, count_column])
        if np.any(part[:, count_column] != count):
            raise ValueError("frequency-dependent sample count")
        part_mean = part[:, mean_column]
        part_m2 = part[:, m2_column]
        if total == 0:
            total = count
            mean = part_mean.copy()
            m2 = part_m2.copy()
            continue
        delta = part_mean - mean
        combined = total + count
        m2 += part_m2 + delta * delta * total * count / combined
        mean += delta * count / combined
        total = combined
    return total, mean, m2


def sample_statistics(values):
    """Return mean, sample standard deviation and SEM over seeds."""
    mean = values.mean(axis=0)
    if len(values) == 1:
        unavailable = np.full(values.shape[1], np.nan)
        return mean, unavailable, unavailable.copy()
    std = values.std(axis=0, ddof=1)
    return mean, std, std / math.sqrt(len(values))


def save_table(path, columns, description, labels):
    """Write one narrow table with explicitly numbered columns."""
    numbered = " ".join(f"{index}.{label}"
                          for index, label in enumerate(labels, 1))
    np.savetxt(path, np.column_stack(columns), fmt="%.17e",
               header=f"{description}\n{numbered}")


def merge_channel(folders, channel, nlat, lfreq, output):
    """Merge spatial populations and preserve seed-to-seed fluctuations."""
    pairs = [read_channel(folder, channel, nlat, lfreq)
             for folder in folders]
    stats_parts = [pair[0] for pair in pairs]
    diag_parts = [pair[1] for pair in pairs]
    omega = stats_parts[0][:, 0]
    for folder, part in zip(folders[1:], stats_parts[1:]):
        if not np.allclose(omega, part[:, 0], rtol=0, atol=1e-13):
            raise ValueError(f"{folder}: incompatible frequency grid")

    count, rho_mean, rho_m2 = merge_population(stats_parts, 2, 3)
    log_count, log_mean, log_m2 = merge_population(stats_parts, 4, 5)
    if count != log_count:
        raise ValueError(f"{channel}: inconsistent linear/log counts")
    spatial_std = np.sqrt(np.maximum(rho_m2 / count, 0.0))
    log_spatial_std = np.sqrt(np.maximum(log_m2 / count, 0.0))

    seed_rho = np.stack([part[:, 2] for part in stats_parts])
    seed_log = np.stack([part[:, 4] for part in stats_parts])
    _, seed_std, seed_sem = sample_statistics(seed_rho)
    _, log_seed_std, log_seed_sem = sample_statistics(seed_log)
    rho_typical = np.exp(log_mean)

    negative_count = sum(part[:, 1] for part in diag_parts)
    floor_count = sum(part[:, 3] for part in diag_parts)
    raw_min = np.minimum.reduce([part[:, 2] for part in diag_parts])
    raw_mean = np.mean([part[:, 4] for part in diag_parts], axis=0)

    save_table(
        output / f"rho_average_{channel}.dat",
        (omega, rho_mean, spatial_std, seed_std, seed_sem),
        f"Arithmetic LDOS statistics for {channel}; {len(folders)} seeds",
        ("omega", "rho_average", "spatial_std", "seed_std", "seed_sem"),
    )
    save_table(
        output / f"rho_typical_{channel}.dat",
        (omega, rho_typical, log_mean, log_spatial_std,
         log_seed_std, log_seed_sem, floor_count / count),
        f"Geometric LDOS statistics for {channel}; floor acts only in log",
        ("omega", "rho_typical", "mean_logrho", "log_spatial_std",
         "log_seed_std", "log_seed_sem", "floor_fraction"),
    )
    save_table(
        output / f"rho_diagnostics_{channel}.dat",
        (omega, negative_count, negative_count / count, raw_min,
         floor_count, floor_count / count, raw_mean),
        f"LDOS numerical diagnostics for {channel}; samples={count}",
        ("omega", "negative_count", "negative_fraction", "raw_min",
         "floor_count", "floor_fraction", "raw_mean"),
    )
    return {
        "channel": channel,
        "omega": omega,
        "mean": rho_mean,
        "typical": rho_typical,
        "spatial_std": spatial_std,
        "seed_sem": seed_sem,
        "negative_fraction": negative_count / count,
        "floor_fraction": floor_count / count,
    }


def plot_channel(result, output, metadata, nseed, dpi, energy_window, plt):
    """Follow the AmorphousSC mean/typical spectrum presentation."""
    omega = result["omega"]
    mean = result["mean"]
    typical = result["typical"]
    std = result["spatial_std"]
    positive = np.concatenate((mean[mean > 0], typical[typical > 0]))
    floor = positive.min() if len(positive) else np.finfo(float).tiny
    lower = np.maximum(mean - std, floor)
    upper = np.maximum(mean + std, floor)

    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    ax.fill_between(omega, lower, upper, color="gray", alpha=0.22,
                    linewidth=0,
                    label=r"$\overline{\rho}\pm\sigma_{\rm spatial}$")
    ax.plot(omega, mean, color="red", linestyle="--", linewidth=1.6,
            label="Arithmetic mean")
    ax.plot(omega, typical, color="black", linewidth=1.5,
            label="Geometric mean")
    ax.set(xlabel=r"$\omega$", ylabel=r"$\rho(\omega)$",
           yscale="log",
           title=(f"{CHANNEL_LABELS[result['channel']]} | "
                  f"{metadata['source']}, {metadata['method']} | "
                  f"{nseed} seeds x {metadata['nlat']} sites | "
                  rf"$\eta={numeric(metadata['ieta']):g}$"))
    if energy_window is not None:
        ax.set_xlim(*energy_window)
    ax.legend()
    ax.grid(alpha=0.18)
    fig.savefig(output / f"spectrum_{result['channel']}.png", dpi=dpi)
    plt.close(fig)


def plot_ratio(results, output, metadata, nseed, dpi, energy_window, plt):
    """Plot the typical-to-average LDOS ratio for all saved channels."""
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for result in results:
        ratio = np.divide(result["typical"], result["mean"],
                          out=np.zeros_like(result["mean"]),
                          where=result["mean"] > 0)
        ax.plot(result["omega"], ratio, linewidth=1.4,
                label=CHANNEL_LABELS[result["channel"]])
    ax.set(xlabel=r"$\omega$",
           ylabel=r"$\rho_{\rm typ}/\rho_{\rm av}$",
           ylim=(0, 1.05),
           title=(f"LDOS heterogeneity | {metadata['source']}, "
                  f"{metadata['method']} | {nseed} seeds"))
    if energy_window is not None:
        ax.set_xlim(*energy_window)
    ax.legend(ncol=2)
    ax.grid(alpha=0.18)
    fig.savefig(output / "rho_typical_ratio.png", dpi=dpi)
    plt.close(fig)


def write_summary(output, metadata, seeds, results, figures):
    """Write compact integral diagnostics and a linked Markdown report."""
    rows = []
    notes = [
        "# Anderson-BHZ spectral ensemble",
        "",
        f"Source: `{metadata['source']}`; method: `{metadata['method']}`; "
        f"seeds: {len(seeds)}; sites per seed: {metadata['nlat']}; "
        f"eta: {numeric(metadata['ieta']):g}.",
        "",
        "Spatial standard deviations describe the pooled local LDOS. "
        "Seed standard deviations and SEM use independent seed means.",
        "",
        "| channel | integral average | integral typical | "
        "max negative fraction | max floor fraction |",
        "|---|---:|---:|---:|---:|",
    ]
    for index, result in enumerate(results, 1):
        average_integral = np.trapezoid(result["mean"], result["omega"])
        typical_integral = np.trapezoid(result["typical"],
                                        result["omega"])
        negative = result["negative_fraction"].max()
        floor = result["floor_fraction"].max()
        rows.append((index, average_integral, typical_integral,
                     negative, floor))
        notes.append(
            f"| {result['channel']} | {average_integral:.8g} | "
            f"{typical_integral:.8g} | {negative:.3e} | {floor:.3e} |")
    save_table(
        output / "spectrum_summary.dat",
        tuple(np.asarray(rows)[:, column] for column in range(5)),
        "Integrated spectral weights and worst numerical diagnostics; "
        "channel IDs follow metadata order",
        ("channel_id", "integral_average", "integral_typical",
         "max_negative_fraction", "max_floor_fraction"),
    )
    if figures:
        notes += ["", "## Figures", ""]
        for result in results:
            name = f"spectrum_{result['channel']}.png"
            notes += [f"### {result['channel']}", "", f"![{name}]({name})", ""]
        notes += ["### Typical-to-average ratio", "",
                  "![rho_typical_ratio.png](rho_typical_ratio.png)", ""]
    (output / "summary.md").write_text("\n".join(notes) + "\n")


def write_ensemble_metadata(path, reference, seeds):
    """Write the common controls and the exact independent seed list."""
    with path.open("w") as stream:
        stream.write("format spectrum_ensemble_v1\n")
        for key in CONSISTENT_KEYS:
            if key != "format":
                stream.write(f"{key} {reference[key]}\n")
        stream.write(f"nseed {len(seeds)}\n")
        stream.write("seeds " + " ".join(map(str, seeds)) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge compact Anderson-BHZ spectrum statistics")
    parser.add_argument("seed_folders", nargs="+", type=Path,
                        help="seed directories containing spectrum_* files")
    parser.add_argument("-o", "--output", type=Path,
                        default=Path("spectrum_ensemble"),
                        help="ensemble output directory")
    parser.add_argument("--aggregate-only", action="store_true",
                        help="write ensemble tables without figures")
    parser.add_argument("--energy-window", nargs=2, type=float,
                        metavar=("EMIN", "EMAX"),
                        help="optional horizontal range for all figures")
    parser.add_argument("--dpi", type=int, default=180,
                        help="figure resolution in dots per inch")
    args = parser.parse_args(argv)
    try:
        if (args.energy_window is not None
                and (not np.isfinite(args.energy_window).all()
                     or args.energy_window[0] >= args.energy_window[1])):
            raise ValueError("invalid --energy-window")
        if args.dpi < 1:
            raise ValueError("--dpi must be positive")
        folders = [folder.resolve() for folder in args.seed_folders]
        reference, seeds = check_metadata(folders)
        nlat = int(reference["nlat"])
        lfreq = int(reference["lfreq"])
        args.output.mkdir(parents=True, exist_ok=True)
        results = [merge_channel(folders, channel, nlat, lfreq,
                                 args.output) for channel in CHANNELS]
        write_ensemble_metadata(args.output / "spectrum_metadata.dat",
                                reference, seeds)
        if not args.aggregate_only:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            for result in results:
                plot_channel(result, args.output, reference, len(seeds),
                             args.dpi, args.energy_window, plt)
            plot_ratio(results, args.output, reference, len(seeds),
                       args.dpi, args.energy_window, plt)
        write_summary(args.output, reference, seeds, results,
                      figures=not args.aggregate_only)
        print(f"Merged {len(seeds)} seeds in {args.output}")
    except (ImportError, OSError, ValueError) as error:
        parser.exit(2, f"merge_spectrum_stats.py: error: {error}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
