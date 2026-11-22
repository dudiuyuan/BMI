#!/usr/bin/env python3
"""DTW with four fixed variants: choose a mode, no parameters to set by hand.

The four variants were chosen from the grid search (tune_dtw.ipynb, no
scaling) on wpj's baseline data and checked on the medium-difficulty data in
zrz/simulation_new; see outputs/final_candidates/. Every mode uses
--scale-floor 0 (no per-gene scaling). "diagonal" means the window that keeps
the warping path near the diagonal.

    python3 dtw_modes.py --mode time --input centered.csv --output outputs/name
    python3 dtw_modes.py --list

The computation is dtw.py itself: this file only fills in the parameters and
records the mode name in run.json (field "variant"). Other dtw.py options can
be passed through unchanged (--clustering, --clustering-seed, --limit, --pair,
--prune ...); the DTW parameters and --scale-floor are fixed by the mode.
"""
import argparse
import json
import sys
from pathlib import Path

import dtw

MODES = {
    'time': {
        'description': '时间惩罚',
        'candidate': 3,
        'params': {'time_penalty': 2},
    },
    'time-dynamics': {
        'description': '时间惩罚 + 表达变化速率',
        'candidate': 8,
        'params': {'extend_penalty': 0.05, 'time_penalty': 0.5, 'dynamics': 0.3},
    },
    'diagonal-time-dynamics': {
        'description': '对角线窗口 + 时间惩罚 + 表达变化速率',
        'candidate': 6,
        'params': {'window': 0.25, 'extend_penalty': 0.05, 'time_penalty': 0.5, 'dynamics': 0.3},
    },
    'diagonal-time': {
        'description': '对角线窗口 + 时间惩罚',
        'candidate': 1,
        'params': {'window': 0.1667, 'extend_penalty': 0.05, 'time_penalty': 2},
    },
}
SCALE_FLOOR = 0.0
FLAGS = {'window': '--window', 'step_penalty': '--step-penalty', 'open_penalty': '--open-penalty',
         'extend_penalty': '--extend-penalty', 'time_penalty': '--time-penalty', 'power': '--power',
         'dynamics': '--dynamics'}
FIXED = set(FLAGS.values()) | {'--scale-floor', '--config-file'}


def config(mode):
    """Full dtw.Config for a mode (anything the mode does not set stays off)."""
    return dtw.Config(**MODES[mode]['params'])


def dtw_arguments(mode):
    """The dtw.py arguments a mode stands for."""
    args = []
    for name, value in MODES[mode]['params'].items():
        args += [FLAGS[name], f'{value:g}']
    return args + ['--scale-floor', f'{SCALE_FLOOR:g}']


def describe():
    lines = []
    for name, m in MODES.items():
        lines.append(f"{name:<24} {m['description']}（候选第 {m['candidate']} 组）")
        lines.append(f"{'':<24} 等同于 python3 dtw.py {' '.join(dtw_arguments(name))}")
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='Modes:\n' + describe())
    parser.add_argument('--mode', choices=list(MODES), help='Which fixed DTW variant to run')
    parser.add_argument('--list', action='store_true', help='Show the modes and their parameters, then exit')
    parser.add_argument('--output', type=Path, help='Result folder (passed to dtw.py)')
    args, rest = parser.parse_known_args(argv)
    if args.list:
        print(describe())
        return
    if not args.mode:
        parser.error('choose --mode (see --list)')
    fixed = [a for a in rest if a.split('=')[0] in FIXED]
    if fixed:
        parser.error(f'{", ".join(fixed)} is fixed by --mode {args.mode}; use dtw.py to set parameters by hand')
    if args.output is None:
        parser.error('--output is required')
    dtw.main(rest + ['--output', str(args.output)] + dtw_arguments(args.mode))
    # Record which variant produced this run (dtw.main has written run.json).
    run_json = args.output/'run.json'
    if run_json.exists():
        meta = json.loads(run_json.read_text(encoding='utf-8'))
        meta['variant'] = {'name': args.mode, 'description': MODES[args.mode]['description'],
                           'candidate': MODES[args.mode]['candidate']}
        run_json.write_text(json.dumps(meta, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
