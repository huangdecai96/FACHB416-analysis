"""Verify day-7 measurements, class-weighted means and final plotting inputs."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data',type=Path,default=Path(__file__).resolve().parent.parent/'FACHB416-source-data')
p.add_argument('--out',type=Path,default=Path('results/cell_dimension_audit'))
a=p.parse_args();D=a.data.resolve();O=a.out.resolve()
if O==D or D in O.parents or O in D.parents:p.error('Output and data folders must be separate')
if O.exists() and any(O.iterdir()):p.error('Output directory must be empty')
cells=pd.read_csv(D/'Morphology/day7_individual_cell_dimensions_with_colony_size.tsv',sep='\t')
assert not cells.duplicated(['sample_id','cell_id']).any()
calc=np.pi*cells.length_um*cells.width_um**2/6
assert np.allclose(calc,cells.volume_um3,rtol=1e-10,atol=1e-9)
classes=pd.read_csv(D/'Morphology/day7_colony_class_volume_means.tsv',sep='\t')
observed=cells.groupby(['sample_id','colony_size_cells']).volume_um3.agg(['count','mean'])
lookup=classes.set_index(['sample_id','colony_size_cells'])
for key,row in lookup.iterrows():
    if row.measured_cell_n>0:
        assert observed.loc[key,'count']==row.measured_cell_n
        assert np.isclose(observed.loc[key,'mean'],row.observed_mean_volume_um3,atol=1e-8,rtol=0)
    elif row.counted_cells>0:
        donor=(row.donor_sample_id,int(row.donor_colony_size_cells))
        assert np.isclose(observed.loc[donor,'mean'],row.used_mean_volume_um3,atol=1e-8,rtol=0)
    else:
        assert row.counted_cells==0 and pd.isna(row.used_mean_volume_um3)
assert np.allclose(classes.counted_cells,classes.counted_particles*classes.colony_size_cells)
weights=classes.counted_cells/classes.groupby('sample_id').counted_cells.transform('sum')
assert np.allclose(weights,classes.cell_fraction,atol=1e-12)
means=(weights*classes.used_mean_volume_um3).groupby(classes.sample_id).sum()
ref=pd.read_csv(D/'Morphology/day7_cell_volume_culture_means.tsv',sep='\t')
ref['sample_id']=ref.condition+ref.replicate.astype(str);ref=ref.set_index('sample_id')
ref['recomputed_class_weighted_um3']=means
ref['difference_um3']=ref.recomputed_class_weighted_um3-ref.cell_volume_um3
assert np.allclose(ref.recomputed_class_weighted_um3,ref.cell_volume_um3,rtol=0,atol=5e-8)
plot=pd.read_csv(D/'Phenotype_expression/plot_replicates.tsv',sep='\t');plot=plot[plot.metric=='day7_cell_volume'].copy()
plot['sample_id']=plot.condition+plot.replicate.astype(str)
assert np.allclose(plot.value,ref.loc[plot.sample_id].cell_volume_um3,atol=1e-8,rtol=0)
O.mkdir(parents=True,exist_ok=True);ref.to_csv(O/'class_weighted_volume_check.tsv',sep='\t')
summary={'measured_cells':len(cells),'cultures':len(ref),'class_weighted_means_verified':True,'plotting_inputs_verified':True,'max_abs_difference_um3':float(abs(ref.difference_um3).max())}
(O/'audit.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
