"""Run with the installed bdmc package, or the project's virtual environment."""
import argparse
from pathlib import Path

from bdmc.experiments import train_calibration

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--model',choices=['brownian','dumbbell','both'],default='both')
parser.add_argument('--without-torch',action='store_true')
args = parser.parse_args()
models = ('brownian','dumbbell') if args.model=='both' else (args.model,)
for model in models:
    report = train_calibration(args.output/model,model=model,include_torch=not args.without_torch)
    print(model,'analytic RMSE',report['analytic_test_rmse'],
          {name:data['metrics']['test']['rmse'] for name,data in report['methods'].items()},flush=True)
