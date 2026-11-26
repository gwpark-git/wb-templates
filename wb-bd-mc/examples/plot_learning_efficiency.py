"""Standalone diagnostic plot from already recorded calibration reports."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--reports',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True,help='fresh filename prefix for SVG and PNG')
args=parser.parse_args()
targets=[args.output.with_suffix(suffix) for suffix in ('.svg','.png')]
if any(path.exists() for path in targets):
    parser.error('plot output exists')
fig,axes=plt.subplots(1,2,figsize=(9,3.7),layout='constrained')
colors={'ridge':'#3366a1','gaussian_process':'#287e58','mlp':'#bb6b22'}
labels={'ridge':'Ridge','gaussian_process':'Gaussian process','mlp':'MLP'}
for ax,model in zip(axes,('brownian','dumbbell')):
    report=json.loads((args.reports/model/'report.json').read_text())
    for method in colors:
        rows=sorted([row for row in report['summary'] if row['method']==method],
                    key=lambda row:row['budget_components'])
        values=[next(iter(row['validation_per_target_rmse'].values())) for row in rows]
        ax.plot([row['budget_components'] for row in rows],[v['mean'] for v in values],
                marker='o',color=colors[method],label=labels[method])
    ax.axhline(report['baseline']['validation_metrics']['rmse'],color='#444444',ls='--',label='Analytic reference')
    ax.set(title=model.capitalize(),xlabel='Independent training conditions',yscale='log',xticks=[1,2,3,5])
    ax.grid(axis='y',which='major',alpha=.2)
    ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('Validation RMSE (normalised squared length)')
axes[1].legend(frameon=False,fontsize=8)
fig.suptitle('Calibration only · 3 seed runs per condition · mean over 3 training subsets',fontsize=10)
for path in targets:
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=180)
plt.close(fig)
print('\n'.join(str(path) for path in targets))
