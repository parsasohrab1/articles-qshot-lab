import argparse
import pathlib

from . import experiment


def main(argv=None):
    ap = argparse.ArgumentParser(prog="qshotlab")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="compare allocation strategies on random Pauli-sum Hamiltonians")
    r.add_argument("--qubits", type=int, default=3)
    r.add_argument("--terms", type=int, default=12)
    r.add_argument("--shots", type=int, nargs="+", default=[500, 2000, 8000])
    r.add_argument("--instances", type=int, default=5)
    r.add_argument("--reps", type=int, default=200)
    r.add_argument("--spread", type=float, default=1.0, help="log-normal spread of the coefficients")
    r.add_argument("--seed", type=int, default=0)
    r.add_argument("--out", default="results")
    a = ap.parse_args(argv)
    d, s = experiment.run(a.qubits, a.terms, a.shots, a.instances, a.reps, a.seed, a.spread)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "instances.csv", index=False)
    s.to_csv(out / "summary.csv", index=False)
    print(s.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
