#!/usr/bin/env python3
"""Track additional generated data against the fixed recorded baseline."""
import argparse
import json
from pathlib import Path
import sys

# Keep this repository tool usable with system Python as well as the venv.
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
from bdmc.storage import allocated_bytes, check, create_data_baseline, check_additional

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('roots',nargs='*',type=Path,help='data roots, only when initializing the baseline')
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--initialize-baseline',action='store_true')
    parser.add_argument('--planned-bytes',type=int,default=0)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.initialize_baseline:
        if not args.roots:
            parser.error('baseline initialization requires explicit data roots')
        create_data_baseline(args.roots,args.baseline)
    elif args.roots:
        parser.error('data roots are recorded in the baseline; do not change them at check time')
    report=check_additional(args.baseline,args.planned_bytes)
    text=json.dumps(report,indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text,end='')
    raise SystemExit(0 if report['status']=='ok' else 1)
