import argparse
from pathlib import Path
from bdmc.experiments import train_network_demo

parser=argparse.ArgumentParser(description='Legacy network compatibility training demo; not scientific validation')
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--executable',type=Path,required=True)
parser.add_argument('--template',type=Path,required=True)
parser.add_argument('--without-torch',action='store_true')
args=parser.parse_args()
report=train_network_demo(args.output,executable=args.executable,template=args.template,
                          include_torch=not args.without_torch)
print('network demo',report['n_runs'],'runs',
      {name:data['metrics']['test']['rmse'] for name,data in report['methods'].items()})
