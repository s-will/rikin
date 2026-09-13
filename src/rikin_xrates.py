#!/usr/bin/env python3
"""Calculate kinetics from weights/pfs of states and transitions by solving
the master equation for Arrhenius rates.

INPUT: file with square matrix of partition functions of 1) transitions and
2) basins (given in diagonal).

OUTPUT: probabilities of states at requested times according to Arrhenius
rates from partition functions; time is increased exponentially.

This is a port of the original Octave implementation (rikin_xrates.m).
"""

import argparse
import struct
import sys

import numpy as np
from scipy.linalg import expm

TOOL_NAME = "XRates"
TOOL_VERSION = "0.2"
TOOL_DESCRIPTION = (
    "Calculate kinetics from weights/pfs of states and transitions by "
    "solving the master equation for Arrhenius rates."
)

# Convergence threshold on the L1 distance between pi(t) and pi(infinity).
CONVERGENCE_THRESHOLD = 5e-3


def tool_header():
    """Reproduce the original's two-line wrapped banner."""
    return "%s %s -- Calculate kinetics from weights/pfs of states and transitions by solving\n  the master equation for Arrhenius rates.\n" % (
        TOOL_NAME,
        TOOL_VERSION,
    )


def read_pfs(filename, binary):
    """Read the square partition-function matrix.

    The binary layout is the one written by rikin_prune --pffile and read by
    the original's `fread(fh,1,"uint64")` / `fread(fh,[dim,dim],"double")`:
    a uint64 dimension header followed by dim*dim native-endian doubles in
    *column-major* (Fortran) order, which is how Octave's fread fills a
    [dim,dim] shape.
    """
    if not binary:
        print("Read text input from file %s" % filename)
        try:
            pfs = np.loadtxt(filename)
        except OSError as e:
            sys.exit("ERROR: cannot read %s: %s" % (filename, e.strerror or e))
        pfs = np.atleast_2d(pfs)
        return pfs, pfs.shape[0]

    print("Read binary input from file %s" % filename)
    try:
        with open(filename, "rb") as fh:
            header = fh.read(8)
            if len(header) < 8:
                sys.exit("ERROR: could not read dimension header from %s."
                         % filename)
            (dim,) = struct.unpack("=Q", header)
            data = np.fromfile(fh, dtype=np.double, count=dim * dim)
    except OSError as e:
        sys.exit("ERROR: cannot read %s: %s" % (filename, e.strerror or e))
    if data.size != dim * dim:
        sys.exit(
            "ERROR: pf file %s holds %d values, expected %d (dim=%d)."
            % (filename, data.size, dim * dim, dim)
        )
    return data.reshape((dim, dim), order="F"), dim


def rates_from_pfs(pfs, basin_pfs):
    """Arrhenius rates R[i,j] = pfs[i,j] / basin_pfs[i], zero diagonal.

    The original writes this as `pfs ./ repmat(basin_pfs,1,dim)`, i.e. a
    row-wise division by the basin weight of the *source* state; the
    broadcast `basin_pfs[:, None]` is the direct equivalent.
    """
    R = pfs / basin_pfs[:, None]
    np.fill_diagonal(R, 0.0)
    return R


def compute_mfpts(R, pi8, dim, absorb, startstate):
    """Mean first passage times, mirroring the original's two branches."""
    print("Compute mean first passage times.")

    # "unit time": one step of the discrete Markov process with transition
    # probability matrix P. The original notes that treekin approximates
    # P ~ I + R, but deliberately uses the exact matrix exponential here.
    timescale = 1.0
    P = expm(timescale * R)

    if absorb > 0:
        a = absorb - 1  # the CLI uses 1-based state indices
        # Q is P without the absorbing state's row and column.
        Q = np.delete(np.delete(P, a, axis=0), a, axis=1)
        # Fundamental matrix N; row sums give absorption times per state.
        N = np.linalg.inv(np.eye(dim - 1) - Q)
        t = N @ np.ones(dim - 1)
        t = np.insert(t, a, 0.0)  # absorption time 0 for the absorbing state
        t = t * timescale
        print(
            "Absorption times (to state %d): %s"
            % (absorb, "".join(" %g" % v for v in t))
        )
    else:
        # Stationary matrix W[i,j] = pi8[j].
        W = np.tile(pi8, (dim, 1))
        # Fundamental matrix Z.
        Z = np.linalg.inv(np.eye(dim) - P + W)
        # m_ij = (z_jj - z_ij) / w_j is the MFPT from i to j.
        M = (np.tile(np.diag(Z), (dim, 1)) - Z) / W
        M = M * timescale
        # Reproduce Octave's disp() under "format short e": 13-character
        # right-aligned columns, exact zeros shown as a bare "0".
        for row in M:
            print("".join(
                "%13s" % ("0" if v == 0 else "%.4e" % v) for v in row
            ))
        print(
            "MFPTs starting from state %d: %s"
            % (startstate, "".join(" %g" % v for v in M[startstate - 1, :]))
        )


def main():
    # The banner is printed before parsing, as in the original.
    print(tool_header())

    p = argparse.ArgumentParser(
        prog="rikin_xrates.py",
        description=TOOL_DESCRIPTION,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("infile", metavar="in",
                   help="Input file. Matrix of transition weights/partition "
                        "functions; the diagonal contains state weights.")
    p.add_argument("--out", required=True,
                   help="Name of output file. The distributions at each time "
                        "point until convergence are written as a table to "
                        "this file.")
    p.add_argument("-v", "--verbose", action="store_true", help="Be verbose.")
    # The original accepts --binary (a no-op, since binary is the default) and
    # selects text input only via --nonbinary; its flag parser also accepted
    # the generic negation --no-binary, kept here as an alias.
    p.add_argument("-b", "--binary", action="store_true",
                   help="Binary input (default; accepted for compatibility).")
    p.add_argument("--nonbinary", "--no-binary", dest="nonbinary",
                   action="store_true", help="Non-binary (text) input.")
    p.add_argument("--t0", type=float, default=1e-1, help="Start time")
    p.add_argument("--t8", type=float, default=1e16, help="End time")
    p.add_argument("--tinc", type=float, default=1.2, help="Time increment")
    p.add_argument("--mode", default="expm", choices=("expm", "diag"),
                   help="Mode (expm or diag)")
    p.add_argument("--p0", type=int, default=2,
                   help="State with initial probability 1 (1-based)")
    p.add_argument("--absorb", type=int, default=-1,
                   help="Absorbing state (1-based)")
    p.add_argument("--mfpts", action="store_true",
                   help="Compute mean first passage times")
    args = p.parse_args()

    verbose = args.verbose
    binary = not args.nonbinary
    starttime, endtime, tinc = args.t0, args.t8, args.tinc
    startstate, absorb, mode = args.p0, args.absorb, args.mode

    if verbose:
        print("Arguments:")
        for k, v in sorted(vars(args).items()):
            print("    %s = %s" % (k, v))

    if verbose:
        print("Load pfs from file %s" % args.infile)

    pfs, dim = read_pfs(args.infile, binary)
    if pfs.shape[0] != pfs.shape[1]:
        sys.exit("ERROR: pf file has to contain a square matrix.")
    if not 1 <= startstate <= dim:
        sys.exit("ERROR: --p0 %d is outside the state range 1..%d."
                 % (startstate, dim))
    if absorb > dim:
        sys.exit("ERROR: --absorb %d is outside the state range 1..%d."
                 % (absorb, dim))

    basin_pfs = np.diag(pfs).copy()

    # ---- symmetry check -------------------------------------------------
    # Rates derived from the original and from the transposed pf matrix are
    # equal iff the input is symmetric; compare via the spectral norm.
    R1 = rates_from_pfs(pfs, basin_pfs)
    R2 = rates_from_pfs(pfs.T, basin_pfs)
    dev_symm = np.linalg.norm(R1 - R2, 2)
    if dev_symm > 0.001:
        print("WARNING: questionable input symmetry of input pfs: "
              "norm_2(R(pfs)-R(t(pfs)))=%g (This value should be almost 0.)"
              % dev_symm)

    # The input has to be symmetric; anything else is wrong input. Force it.
    pfs = (pfs + pfs.T) / 2.0

    # ---- rate matrix ----------------------------------------------------
    R = rates_from_pfs(pfs, basin_pfs)
    # Diagonal is -rowsum, so that each row of R sums to zero.
    np.fill_diagonal(R, -R.sum(axis=1))

    # ---- initial and stationary distribution ----------------------------
    pi0 = np.zeros(dim)
    pi0[startstate - 1] = 1.0
    if verbose:
        print("PI_0: %s" % "".join(" %g" % v for v in pi0))

    pi8 = basin_pfs / basin_pfs.sum()  # pi after infinite time

    if absorb > 0:
        a = absorb - 1
        pi8 = np.zeros(dim)
        pi8[a] = 1.0
        R[a, :] = 0.0  # no outflow from the absorbing state

    if verbose:
        print("PI_8: %s" % "".join(" %g" % v for v in pi8))
        print()
        print("Compute distributions at times %e..%e (until convergence)"
              % (starttime, endtime))

    if args.mfpts:
        compute_mfpts(R, pi8, dim, absorb, startstate)

    # From here on we need R[i,j] = rate from j to i.
    R = R.T

    # ---- optional precomputation for diagonalization mode ---------------
    if mode == "diag":
        # Symmetrize by the similarity transform diag(sqrt(pi8))^-1 R
        # diag(sqrt(pi8)), valid for rate matrices in detailed balance.
        # ATTENTION: this requires the correct pi8.
        sqrt_pi8 = np.sqrt(pi8)
        # Translate by +I first, so the transformed matrix stays positive
        # definite for the eigensolver, then undo it on the eigenvalues.
        symmR = np.eye(dim) + R
        symmR = (symmR / sqrt_pi8[:, None]) * sqrt_pi8[None, :]

        # eigh (not eig) is the correct solver for a symmetric matrix: it
        # returns real eigenvalues and a genuinely orthonormal eigenvector
        # basis, which is what makes the transpose below a valid inverse.
        eigvals, eigvecs = np.linalg.eigh(symmR)
        eigvecs_inv = eigvecs.T
        eigvals = eigvals - 1.0  # compensate for the +I translation

        # Control: the code below uses eigvecs.T as the inverse of eigvecs,
        # which is only valid if the eigenvectors are orthonormal. This holds
        # for the symmetrized (detailed-balance) rate matrix, but can degrade
        # numerically for ill-conditioned input -- in which case every
        # distribution computed from the decomposition is wrong. Measuring how
        # far eigvecs.T @ eigvecs is from the identity tests exactly that
        # assumption. Unlike a residual involving the rates themselves, this is
        # scale-free: it is near 0 whenever the decomposition is sound,
        # regardless of how large the rates are.
        dev_diag = np.linalg.norm(eigvecs_inv @ eigvecs - np.eye(dim), 2)
        print("Control diagonalization: %g (value should be almost 0)"
              % dev_diag)

        pre_left = sqrt_pi8[:, None] * eigvecs
        pre_right = eigvecs_inv @ (pi0 / sqrt_pi8)

    # ---- main loop ------------------------------------------------------
    step = 0
    time = starttime
    with open(args.out, "w") as fout:
        while time < endtime:
            if mode == "diag":
                pi = pre_left @ (np.exp(time * eigvals) * pre_right)
            else:
                pi = expm(time * R) @ pi0

            fout.write("%e" % time)
            fout.write("".join(" %g" % v for v in pi))
            fout.write("\n")

            # Converged once pi(t) is close to pi(infinity) in L1.
            if np.abs(pi - pi8).sum() < CONVERGENCE_THRESHOLD:
                print("Convergence at time %g." % time)
                break

            time *= tinc
            step += 1

    # Note: as in the original, the step counter is not incremented for the
    # iteration that detects convergence, so this is one less than the
    # number of lines written.
    print("Computed distributions at %d time points" % step)


if __name__ == "__main__":
    main()
