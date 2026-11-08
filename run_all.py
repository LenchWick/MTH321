"""Run Group D's Project 2 tests, experiments, figure export and verification."""
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true', help='check existing outputs without remeasuring experiments')
    parser.add_argument('--output-dir', type=Path, help='write a separate rerun without replacing the published reference outputs')
    args = parser.parse_args()
    out = args.output_dir.resolve() if args.output_dir else ROOT/'code/project2/outputs'
    figures = out/'figures' if args.output_dir else ROOT/'figures/project2'
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['MTH321_OUTPUT_DIR'] = str(out)
    env.setdefault('MPLCONFIGDIR', str(ROOT/'.cache/matplotlib'))
    env.setdefault('XDG_CACHE_HOME', str(ROOT/'.cache'))
    tests = subprocess.run([sys.executable, str(ROOT/'code/project2/test_solver.py')], cwd=ROOT, env=env,
                           text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (out/'unit_tests.txt').write_text(tests.stdout)
    print(tests.stdout, end='')
    tests.check_returncode()
    if not args.check_only:
        subprocess.run([sys.executable, str(ROOT/'code/project2/run_all.py')], cwd=ROOT, env=env, check=True)
        figures.mkdir(parents=True, exist_ok=True)
        for item in out.iterdir():
            if item.suffix.lower() in {'.png', '.pdf', '.gif'}:
                shutil.copy2(item, figures/item.name)
    subprocess.run([sys.executable, str(ROOT/'code/project2/verify_results.py'),
                    '--output-dir', str(out), '--figure-dir', str(figures)], cwd=ROOT, env=env, check=True)
    print(f'Project 2 verified. Data: {out}\nFigures: {figures}')


if __name__ == '__main__':
    main()
