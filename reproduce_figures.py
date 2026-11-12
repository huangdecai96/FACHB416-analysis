"""Render selected final figure definitions into separate new output directories."""
from pathlib import Path
import argparse, subprocess
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--rscript',default='Rscript')
p.add_argument('--figures',nargs='+',default=['Figure4','Figure5'])
a=p.parse_args();repo=Path(__file__).resolve().parent
allowed=['Figure1','Figure2','Figure3','Figure4','Figure5','FigureS2','FigureS5','FigureS6','FigureS7','FigureS8','FigureS9','FigureS10','FigureS11']
if not all(x in allowed for x in a.figures):p.error('Unsupported figure; see figures/Non_R_figures.md')
if a.out.exists() and any(a.out.iterdir()):p.error('Use a new empty output directory')
a.out.mkdir(parents=True,exist_ok=True)
for f in a.figures:
    subprocess.run([a.rscript,str(repo/'figures'/f'{f}.R'),str(a.data.resolve()),str((a.out/f).resolve())],check=True)
print('Completed:',', '.join(a.figures))
