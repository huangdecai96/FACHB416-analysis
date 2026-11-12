"""Replay final plotted statistics from the source-data archive (Python 3; numpy, pandas, scipy).

Run: python reproduce_statistics.py --data ../FACHB416-source-data
Exhaustive six-condition correlations use all 720 permutations. Pairwise Welch
tests use one BH family of 24; interactions one family of 4; correlations 28.
This script does not refit DESeq2 or the upstream genome assembly.
"""
from pathlib import Path
from itertools import permutations
import json
import numpy as np
import pandas as pd
from scipy import stats

import argparse
parser=argparse.ArgumentParser(description="Replay the final FACHB-416 numerical analyses")
parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parent.parent/"FACHB416-source-data")
parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent/"results/statistics")
args=parser.parse_args(); D=args.data.resolve(); O=args.out.resolve()
if not (D/"figure_data_map.tsv").is_file(): parser.error("--data must point to the extracted FACHB416-source-data folder")
if O==D or D in O.parents or O in D.parents: parser.error("Output and source-data folders must be separate")
if O.exists() and any(O.iterdir()): parser.error("Output folder is nonempty; select a new --out folder")
O.mkdir(parents=True,exist_ok=True)
rd=lambda p:pd.read_csv(D/p,sep='\t')
def bh(p):
    p=np.asarray(p,dtype=float);order=np.argsort(p);ans=np.empty(len(p));ans[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return ans
def save(df,name):df.to_csv(O/name,sep='\t',index=False,float_format='%.15g')
checks={}
def close(name,actual,expected,atol=1e-10):
    a=np.asarray(actual,dtype=float);b=np.asarray(expected,dtype=float)
    error=float(np.max(np.abs(a-b)));checks[name]={'max_abs_error':error,'passed':bool(np.allclose(a,b,rtol=0,atol=atol))}
    assert checks[name]['passed'],(name,error)

m=rd('Expression/module_sample_scores.tsv');membership=rd('Expression/gene_set_membership.tsv')
vst=pd.read_csv(D/'RNA/paper2_finalref_vst_matrix.csv',index_col=0)
replay=[]
for module,group in membership.groupby('score_module'):
    genes=group.loc[group.detected==1,'gene_id'];assert genes.isin(vst.index).all()
    scores=vst.loc[genes].mean(axis=0)
    observed=m[m.score_module==module].set_index('sample').score
    close('VST_mean_'+module,scores.loc[observed.index],observed)
    replay.extend(dict(sample=s,score_module=module,score=x) for s,x in scores.items())
save(pd.DataFrame(replay),'module_scores_from_VST.tsv')

contrasts={'HN_vs_LN':('HN','LN'),'HU_vs_LU':('HU','LU'),'HA_vs_LA':('HA','LA'),'HU_vs_HN':('HU','HN'),'HA_vs_HN':('HA','HN'),'HA_vs_HU':('HA','HU')}
rows=[];interactions=[]
for module,d in m.groupby('score_module'):
    assert len(d)==18 and (d.groupby('condition').size()==3).all()
    for contrast,(a,b) in contrasts.items():
        x=d.loc[d.condition==a,'score'];y=d.loc[d.condition==b,'score'];t=stats.ttest_ind(x,y,equal_var=False)
        rows.append(dict(score_module=module,contrast=contrast,delta=x.mean()-y.mean(),t=t.statistic,df=t.df,p=t.pvalue))
    # Compare additive and full source-by-level linear models, F(2,12).
    source=d.condition.str[1];level=(d.condition.str[0]=='H').astype(float).to_numpy()
    u=(source=='U').astype(float).to_numpy();a=(source=='A').astype(float).to_numpy()
    X0=np.column_stack([np.ones(18),u,a,level]);X1=np.column_stack([X0,u*level,a*level]);y=d.score.to_numpy()
    sse=lambda X:np.sum((y-X@np.linalg.lstsq(X,y,rcond=None)[0])**2)
    f=((sse(X0)-sse(X1))/2)/(sse(X1)/12)
    interactions.append(dict(score_module=module,F_interaction=f,p_interaction=stats.f.sf(f,2,12)))
out=pd.DataFrame(rows);out['q_BH24']=bh(out.p);out['stars']=np.select([out.q_BH24<.001,out.q_BH24<.01,out.q_BH24<.05],['***','**','*'],default='ns')
ref=rd('Expression/Pairwise_tests_BH24.tsv');cmp=out.merge(ref,on=['score_module','contrast'],suffixes=('_new','_ref'))
for col in ['p','q_BH24','delta']:close('Welch_'+col,cmp[col+'_new'],cmp[col+'_ref'])
save(out,'Pairwise_tests_BH24.tsv')
it=pd.DataFrame(interactions);it['q_interaction_BH4']=bh(it.p_interaction);cmp=it.merge(rd('Expression/Interaction_BH4_verified.tsv'),on='score_module',suffixes=('_new','_ref'))
close('Interaction_BH4',cmp.q_interaction_BH4_new,cmp.q_interaction_BH4_ref);save(it,'Interaction_BH4.tsv')

c=rd('Phenotype_expression/figure5_correlation_input.tsv');out=[];perm=np.array(list(permutations(range(6))))
for (module,metric),d in c.groupby(['score_module','metric']):
    assert len(d)==6 and d.condition.nunique()==6
    x=stats.rankdata(d['mean_score']);y=stats.rankdata(d['mean']);x=x-x.mean();y=y-y.mean()
    den=np.linalg.norm(x)*np.linalg.norm(y);rho=x@y/den
    null=y[perm]@x/den;p=float(np.mean(abs(null)>=abs(rho)-1e-12))
    out.append(dict(score_module=module,metric=metric,rho=rho,p_exact=p,permutations=720))
out=pd.DataFrame(out);out['q_BH28']=bh(out.p_exact);cmp=out.merge(rd('Phenotype_expression/figure5_correlations.tsv'),on=['score_module','metric'],suffixes=('_new','_ref'))
for col in ['rho','p_exact','q_BH28']:close('Spearman_'+col,cmp[col+'_new'],cmp[col+'_ref'])
save(out,'Correlations_BH28.tsv')

phen=rd('Phenotype_expression/plot_replicates.tsv');anova=[]
for metric,d in phen.groupby('metric'):
    t=stats.f_oneway(*[g.value for _,g in d.groupby('condition')]);anova.append(dict(metric=metric,anova_p=t.pvalue))
anova=pd.DataFrame(anova);cmp=anova.merge(rd('Phenotype_expression/figure4_global_tests.tsv'),on='metric',suffixes=('_new','_ref'));close('Phenotype_ANOVA',cmp.anova_p_new,cmp.anova_p_ref);save(anova,'Phenotype_ANOVA.tsv')
raw=rd('Morphology/colony_counts_all_days.tsv');assert len(raw)==162 and (raw.groupby(['condition','day']).size()==3).all()
cpp=(raw.n1+2*raw.n2+4*raw.n4+8*raw.n8)/(raw.n1+raw.n2+raw.n4+raw.n8);close('CPP_recalculation',cpp,raw.CPP)
day7=raw[raw.day==7].merge(phen[phen.metric=='day7_cpp'],on=['condition','replicate']);close('CPP_final_plot',day7.CPP,day7.value)
close('Fraction_sum',raw.filter(like='_cell_fraction').sum(axis=1),np.ones(len(raw)))
summary=raw.groupby(['condition','day']).CPP.agg(['count','mean','std']).reset_index();summary['se']=summary['std']/np.sqrt(summary['count']);save(summary,'CPP_summary.tsv')

degs=[]
for con in contrasts:
    d=pd.read_csv(D/f'RNA/paper2_finalref_condition_{con}.csv');s=d[(d.padj<.05)&(abs(d.log2FoldChange)>=1)]
    degs.append(dict(contrast=con,tested_genes=len(d),sig_genes=len(s),up_genes=int((s.log2FoldChange>0).sum()),down_genes=int((s.log2FoldChange<0).sum())))
degs=pd.DataFrame(degs);ref=pd.read_csv(D/'RNA/paper2_finalref_condition_DEG_summary_padj0.05_lfc1.csv');cmp=degs.merge(ref,on='contrast',suffixes=('_new','_ref'))
for col in ['tested_genes','sig_genes','up_genes','down_genes']:close('DEG_'+col,cmp[col+'_new'],cmp[col+'_ref'],atol=0)
save(degs,'DEG_summary.tsv')

# Re-aggregate the final family-level six-model consensus for every branch.
fam=rd('Comparative_genomics/replicate_stable_six_model_family_consensus.tsv')
stable=fam[fam.replicate_stable_all_six.astype(str).str.lower()=='true']
agg=stable.groupby(['dataset','branch_label','dominant_direction']).size().unstack(fill_value=0).reset_index()
agg=agg.rename(columns={'contraction':'stable_six_contracted','expansion':'stable_six_expanded'})
branch=rd('Comparative_genomics/replicate_stable_terminal_branch_summary.tsv')
cmp=branch.merge(agg,on=['dataset','branch_label'],how='left',suffixes=('_ref','_new')).fillna(0)
for col in ['stable_six_contracted','stable_six_expanded']:close('CAFE_'+col,cmp[col+'_new'],cmp[col+'_ref'],atol=0)
save(agg,'CAFE_consensus_summary.tsv')
checks['counts']={'CPP':len(raw),'RNA_libraries':len(vst.columns),'Welch_tests':24,'interaction_tests':4,'correlations':28,'all_checks_passed':True}
(O/'verification.json').write_text(json.dumps(checks,indent=2),encoding='utf8');print(json.dumps(checks,indent=2))
