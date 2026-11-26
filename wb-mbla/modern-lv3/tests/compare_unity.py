"""Compare the injected unity law with the Level 2 built-in unity closures."""

from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def run_case(executable, directory, case):
    """Run a controlled fixture in its own directory and retain all three tables."""
    directory.mkdir()
    (directory / "case.inp").write_text(case, encoding="utf-8")
    result = subprocess.run([executable, "case.inp"], cwd=directory,
                            capture_output=True, text=True, timeout=20)
    if result.returncode or "MBLA_STATUS=converged" not in result.stdout:
        raise RuntimeError(result.stdout + result.stderr)
    return {name: (directory / name).read_bytes() for name in
            ("result.dat", "result.dat.Fval", "result.dat.Rout")}


def main():
    level2, level3, fixture = sys.argv[1:]
    live = Path(fixture).read_text(encoding="utf-8")
    builtin = live.replace("fcn_Dc_div_D0=LIVE", "fcn_Dc_div_D0=CT").replace(
        "fcn_eta_div_eta0=LIVE", "fcn_eta_div_eta0=CT")
    with TemporaryDirectory(prefix="mbla-lv3-unity-") as directory:
        root = Path(directory)
        expected = run_case(level2, root / "builtin", builtin)
        observed = run_case(level3, root / "provider", live)
    for name, data in expected.items():
        if not data or observed[name] != data:
            raise RuntimeError(f"Unity transport output differs for {name}")
        print(f"{name}: {len(data.splitlines())} lines, exact byte equality")


if __name__ == "__main__":
    main()
