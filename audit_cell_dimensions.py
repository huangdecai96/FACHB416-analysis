"""Verify measured dimensions and compare direct means with archived plot inputs."""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data',type=Path,default=Path(__file__).resolve().parent.parent/'FACHB416-source-data')
p.add_argument('--out',type=Path,default=Path('results/cell_dimension_audit'))
a=p.parse_args()
if a.out.exists() and any(a.out.iterdir()):p.error('Output directory must be empty')
a.out.mkdir(parents=True,exist_ok=True)
cells=pd.read_csv(a.data/'Morphology/day7_individual_cell_dimensions.tsv',sep='\t')
counts=pd.read_csv(a.data/'Morphology/day7_measurement_counts.tsv',sep='\t')
assert not cells.duplicated(['sample_id','cell_id']).any()
calc=np.pi*cells.length_um*cells.width_um**2/6
assert np.allclose(calc,cells.volume_um3,rtol=1e-10,atol=1e-9)
actual=cells.groupby('sample_id').size()
expected=counts.set_index('sample_id').n_cells
assert actual.sort_index().equals(expected.sort_index())
means=cells.assign(recalculated_volume_um3=calc).groupby('sample_id').recalculated_volume_um3.mean()
reference=pd.read_csv(a.data/'Morphology/day7_cell_volume_culture_means.tsv',sep='\t')
reference['sample_id']=reference.condition+reference.replicate.astype(str)
reference=reference.set_index('sample_id')
reference['direct_cell_mean_um3']=means
reference['direct_minus_archived_um3']=reference.direct_cell_mean_um3-reference.cell_volume_um3
reference.to_csv(a.out/'direct_vs_archived.tsv',sep='\t')
agrees=bool(np.allclose(reference.direct_cell_mean_um3,reference.cell_volume_um3,rtol=1e-10,atol=1e-8))
summary=dict(cell_count=len(cells),culture_count=len(actual),counts_verified=True,
    ellipsoid_formula_verified=True,archived_means_match_direct_means=agrees,
    manuscript_values_replaced=False,
    note='Direct means and archived plotting values are compared without changing either input.')
(a.out/'audit.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
print(json.dumps(summary,indent=2))
