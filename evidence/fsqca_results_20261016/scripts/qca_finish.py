from qca_setup import *
import pandas as pd,numpy as np,itertools,platform,sys
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
stage=Path(__file__).parent
d=pd.read_csv(OUT/'fsqca_calibrated_dataset.csv');meta=pd.read_csv(EV/'scenario_metadata.csv');tt=pd.read_csv(EV/'Python_truth_tables.csv',dtype={'corner':str}).fillna({'countries':'','contradictory_cases':''});nec=pd.read_csv(EV/'Python_necessity.csv');sol=pd.read_csv(EV/'Python_solutions.csv').fillna({'expression':'','counterfactual_corners':''});paths=pd.read_csv(EV/'Python_pathways.csv');support=pd.read_csv(EV/'primary_pathway_case_support.csv')
plan=json.loads((EV/'approved_execution_plan_frozen.json').read_text());freeze=json.loads((EV/'calibration_frozen_receipt.json').read_text())
def primary(x):return x[x.scenario.str.startswith('primary')].copy()
checks=[]
def check(label,error,tol=1e-10):checks.append(dict(check=label,max_absolute_difference=float(error),tolerance=tol,passed=bool(error<=tol)))
rdisplay=pd.read_csv(EV/'R_truth_tables.csv',dtype={'corner':str,'OUT':str});rtt=pd.read_csv(EV/'R_truth_tables_unmasked.csv',dtype={'corner':str}).merge(rdisplay[['scenario','corner','OUT']],on=['scenario','corner']);rn=pd.read_csv(EV/'R_necessity.csv');rs=pd.read_csv(EV/'R_solutions.csv').fillna({'expression':''});rp=pd.read_csv(EV/'R_pathways.csv')
j=tt.merge(rtt,on=['scenario','corner'],suffixes=('_Python','_R'));assert len(j)==len(tt)
for p,r in [('n','n'),('consistency','incl'),('PRI','PRI')]:check('All truth-table '+p,np.max(abs(j[p+'_Python' if p in rtt.columns else p]-j[r+'_R' if r in tt.columns else r])))
check('All truth-table classifications',np.sum(j.OUT_Python!=j.OUT_R),0)
j=nec.merge(rn,on=['scenario','condition'],suffixes=('_Python','_R'));assert len(j)==len(nec)
for c in ['consistency','coverage','RoN']:check('All necessity '+c,np.max(abs(j[c+'_Python']-j[c+'_R'])))
def canonical(x):return ' + '.join(sorted('*'.join(sorted(t.strip().split('*'))) for t in x.split('+'))) if x else ''
sol['canonical']=sol.expression.map(canonical);rs['canonical']=rs.expression.map(canonical)
# R records one 'none' row; independent Python records both unavailable solution types.
js=sol[sol.expression!=''].merge(rs[rs.expression!=''],on=['scenario','type','canonical'],suffixes=('_Python','_R'))
check('Number of independently matched nonempty solutions',abs(len(js)-len(sol[sol.expression!=''])),0)
for c in ['consistency','PRI','coverage']:check('All solution '+c,np.max(abs(js[c+'_Python']-js[c+'_R'])))
jp=paths.merge(rp,on=['scenario','type','model','pathway'],suffixes=('_Python','_R'))
check('Number of independently matched pathways',abs(len(jp)-len(paths)),0)
for c in ['consistency','PRI','raw_coverage','unique_coverage']:check('All pathway '+c,np.max(abs(jp[c+'_Python']-jp[c+'_R'])))
check('R minimization failures',0 if not (EV/'R_errors.csv').exists() else len(pd.read_csv(EV/'R_errors.csv')),0)
check('Calibrated file unchanged since pre-truth-table freeze',int(sha(OUT/'fsqca_calibrated_dataset.csv')!=freeze['calibrated_sha256']),0)
validation=pd.DataFrame(checks);validation.to_csv(EV/'independent_validation.csv',index=False)
if not validation.passed.all():print(validation.to_string(index=False));raise AssertionError('Independent verification failed: resolve before reporting')
pt=primary(tt);pn=primary(nec);ps=primary(sol);pp=primary(paths)
aug=meta.merge(tt.groupby('scenario').agg(positive_rows=('OUT',lambda x:sum(x=='1')),frequency_remainders=('OUT',lambda x:sum(x=='?')),max_PRI=('PRI','max'),max_consistency=('consistency','max')).reset_index(),on='scenario')
complex_sol=sol[sol.type=='complex'][['scenario','expression','consistency','PRI','coverage','union_pass']].rename(columns={c:'solution_'+c for c in ['expression','consistency','PRI','coverage','union_pass']})
aug=aug.merge(complex_sol,on='scenario');aug['primary_expression']=aug.outcome.map({'Ypc':'','Yabs':'FD + RQ'});aug['same_expression_as_primary']=aug.solution_expression==aug.primary_expression
aug.to_csv(EV/'sensitivity_summary.csv',index=False)
coverage=pd.read_csv(ROOT/'fsqca_preparation_20261016/audit_evidence/case_coverage_by_income.csv');coverage['case_retention_share']=coverage.eligible_cases/coverage.economies
allcases=pd.read_csv(ROOT/'fsqca_preparation_20261016/audit_evidence/all_economy_case_eligibility.csv')
groups=['Low income','Lower middle income','Upper middle income','High income']
corner=[''.join(map(str,b)) for b in itertools.product([0,1],repeat=3)]
d['corner']=d[['FD','RQ','GDP']].gt(.5).astype(int).astype(str).agg(''.join,axis=1)
grp=[];cg=[]
for g in groups:
 z=d[d.income_group_2018==g]
 grp.append(dict(income_group=g,cases=len(z),high_Ypc=int((z.Ypc>.5).sum()),high_Yabs=int((z.Yabs>.5).sum()),observed_configurations=z.corner.nunique(),stable_income_2018_2023=int(z.stable_income_2018_2023.sum()),country_codes=';'.join(z.country_code),country_names=';'.join(z.country_name)))
 for c in corner:
  zz=z[z.corner==c];cg.append(dict(income_group=g,corner=c,cases=len(zz),high_Ypc=int((zz.Ypc>.5).sum()),high_Yabs=int((zz.Yabs>.5).sum()),country_codes=';'.join(zz.country_code)))
grp=pd.DataFrame(grp);cg=pd.DataFrame(cg)
pg=support.groupby(['scenario','type','model','pathway','income_group']).agg(cases=('country_code','size'),pathway_owners=('pathway_owner','sum'),supporting_high_outcome=('supports_high_outcome','sum'),contradictions=('contradiction','sum'),membership_mass=('pathway_membership','sum')).reset_index()
near=[];anchors=[]
bands={'FD':5,'RQ':.1,'GDP':800,'Ypc':1,'Yabs':2.5}
for s,(raw,a) in plan['sets'].items():
 anchors.append(dict(set=s,raw_variable=raw,full_out=a[0],crossover=a[1],full_in=a[2],full_out_membership=.05,crossover_membership=.5,full_in_membership=.95,near_band=bands[s],exact_crossovers=int((d[s]==.5).sum()),members_above_crossover=int((d[s]>.5).sum()),minimum_membership=d[s].min(),maximum_membership=d[s].max()))
 for _,r in d[abs(d[raw]-a[1])<=bands[s]].iterrows():near.append(dict(set=s,country_code=r.country_code,country_name=r.country_name,income_group=r.income_group_2018,raw_value=r[raw],crossover=a[1],distance=r[raw]-a[1],membership=r[s]))
anchors=pd.DataFrame(anchors);near=pd.DataFrame(near)
sources=pd.DataFrame([
 dict(source='Approved country cases',definition='One country per case; complete 2018–2020 conditions and 2021–2023 financing/population; original raw values and flags retained.',location=str(SOURCE),sha256=sha(SOURCE)),
 dict(source='IRENA Renewable Energy Statistics 2025',definition='Frozen PDF-validated selected public renewable-energy financing commitments, not disbursements, private investment or all public finance. Constant 2022 USD; raw national annual amounts in millions.',location='https://www.irena.org/Publications/2025/Jul/Renewable-energy-statistics-2025',sha256='See preserved project financing source audit; no new financing vintage used'),
 dict(source='WDI bank private credit',definition='Domestic credit to private sector by banks, % GDP (FD.AST.PRVT.GD.ZS); depth, not overall financial development or stability. Observation-status not disclosed is not proof of unestimated observations.',location='https://data.worldbank.org/indicator/FD.AST.PRVT.GD.ZS',sha256='Existing frozen source provenance retained in case file'),
 dict(source='WGI Regulatory Quality',definition='Governance estimate; uncertainty and source revisions limit sharp calibration interpretation. Zero is reference-scale location, not proven adequacy boundary.',location='https://www.worldbank.org/en/publication/worldwide-governance-indicators',sha256='Existing frozen source provenance retained'),
 dict(source='QCA 3.25 R implementation',definition='Primary truthTable/minimize/pof on exact custom R logistic memberships; no package tail clipping. Independent NumPy sums and exhaustive Boolean cube minimizer.',location='https://cran.r-project.org/package=QCA',sha256=sha(RAW/'QCA_3.25.zip')),
 dict(source='QCA method documentation',definition='Calibration, truth-table, fit and minimization references; author source downloads are inspectable and may differ from released binaries.',location='https://bookdown.org/dusadrian/QCAbook/',sha256='Installed binary results independently verified')])
discrep=pd.read_csv(EV/'R_calibration_checks.csv');discrep['explanation']='QCA::calibrate clips extreme scores below .0001 / above .9999 to 0/1. Approved custom R plogis preserves unrounded logistic scores, independently matched in Python; package clip discrepancy retained.'
rules=pd.DataFrame([dict(rule=k,definition=('UTC '+str(v) if k.endswith('_utc') else str(v))) for k,v in plan.items() if k not in ['sets','input_sha256','sensitivity_operationalization']]);sens_rules=pd.DataFrame([dict(rule=k,definition=str(v)) for k,v in plan['sensitivity_operationalization'].items()])
def clean_rows(df):
 df=df.replace([np.inf,-np.inf],np.nan).astype(object);df=df.where(pd.notna(df),None);return df.to_dict(orient='records')
def sheet(name,df,note):
 if 'corner' in df.columns:
  df=df.copy();df['corner']='C'+df['corner'].astype(str);note+=' Corner C000 = FD/RQ/GDP: 0/0/0.'
 return dict(name=name,note=note,columns=list(df.columns),rows=clean_rows(df))
primarycols=['country_code','country_name','income_group_2018','FD','RQ','GDP','Ypc','Yabs','corner']
rawcols=['country_code','country_name','income_group_2018']+[v[0] for v in plan['sets'].values()]+['financing_ratio_of_sums_per_capita','publisher_estimation_status','source_vintage']
outcomes=pd.DataFrame([dict(outcome=y,concept=('Selected-public commitment intensity per resident' if y=='Ypc' else 'Absolute selected-public commitment volume'),cases=111,high_membership_cases=int((d[y]>.5).sum()),positive_configurations=int((pt[pt.scenario=='primary_'+y].OUT=='1').sum()),sufficient_solutions=int((ps[(ps.scenario=='primary_'+y)&(ps.type=='complex')].expression!='').sum()),necessity_screen_passes=int(pn[pn.scenario=='primary_'+y].screen_pass.sum())) for y in ['Ypc','Yabs']])
payload={
 'fsqca_membership_and_calibration.xlsx':[sheet('Calibrated cases',d[primarycols],'Frozen exact logistic memberships; one country/case; all historical 2018 income groups.'),sheet('Approved anchors',anchors,'Raw-unit anchors map to .05/.50/.95; no epsilon shifts or exact 0/1 clipping.'),sheet('Raw case means',d[rawcols],'Complete three-year means from unchanged approved case file; source and estimation flags retained.'),sheet('Near crossover cases',near,'Review bands fixed in original protocol; no outcome-driven reassignment.'),sheet('Independent checks',validation,'All R–Python comparisons; tolerance 1e-10; classifications and identifiers exact.'),sheet('Package calibration difference',discrep,'Explicit difference: package tail clipping versus approved exact piecewise logistic.'),sheet('Definitions and sources',sources,'Institutional data are inherited frozen sources; no financing vintages spliced.')],
 'fsqca_necessity_results.xlsx':[sheet('Primary necessity',pn,'All conditions and negations; screen .90; no claim of causal necessity.'),sheet('Necessity sensitivity',nec.merge(meta[['scenario','outcome','kind','necessity']],on='scenario'),'Complete local criterion/calibration/deletion results, including failed claims.')],
 'fsqca_truth_tables.xlsx':[sheet('Primary truth tables',pt,'All eight corners per outcome; n>=3, consistency>=.80 and PRI>=.65; frequency remainders labelled ?.'),sheet('Truth table case assignment',d[primarycols],'Exact .5 condition scores belong to no corner; none occur in primary data.'),sheet('Rules',rules,'Thresholds and counterfactual rules frozen before truth tables; intermediate withheld.')],
 'fsqca_sufficiency_solutions.xlsx':[sheet('Primary solution overview',outcomes,'Per-capita and absolute outcomes are distinct. No per-capita sufficient configuration.'),sheet('Primary solutions',ps,'All competing minimum-term covers; union fit assessed separately; no intermediate solution or core/peripheral attribution.'),sheet('Primary pathways',pp,'Raw and unique coverage; codes for owners, high-outcome support and contradictions.'),sheet('Pathway case evidence',support,'Every primary pathway/country membership; zero ownership is not exclusion from fuzzy fit sums.')],
 'fsqca_sensitivity_analyses.xlsx':[sheet('Sensitivity summary',aug,'All 275 frozen scenarios including 222 country deletions; no preferred calibration selected.'),sheet('All truth tables',tt,'Every configuration for every scenario; failed configurations and remainders retained.'),sheet('All solution models',sol,'Complex and parsimonious; empty solutions retained, not relabelled as positive evidence.'),sheet('All pathway fits',paths,'All solution terms, coverage and actual case support.'),sheet('Prespecified local checks',sens_rules,'One factor at a time; local numeric necessity/consistency/PRI alternatives operationalized before calibration.')],
 'fsqca_income_group_comparison.xlsx':[sheet('Group support',grp,'2018 contextual groups, not extra conditions; no separate confirmatory subgroup minimizations.'),sheet('Source selection',coverage,'217-economy reference frame includes territories; high income only 11/80 retained.'),sheet('Configuration by group',cg,'All eight configurations within each income group, including zero observed cells.'),sheet('Pathway by group',pg,'Ownership, supporting cases, contradictions and fuzzy membership mass; descriptive only.'),sheet('High income cases',d[d.income_group_2018=='High income'][primarycols+['stable_income_2018_2023']],'Selected high-income recipients only; cannot generalize to the high-income universe.'),sheet('All economy eligibility',allcases[['country_code','country_name','entity_type','income_group_2018','eligible','exclusion_reasons']],'Original 217-economy eligibility audit; no missing observations recoded to zero.')]
}
(EV/'workbook_payload.json').write_text(json.dumps(payload,allow_nan=False),encoding='utf-8')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','figure.dpi':140})
def save(fig,name):
 fig.savefig(FIG/(name+'.png'),dpi=200);fig.savefig(FIG/(name+'.pdf'));plt.close(fig)
fig,axs=plt.subplots(2,3,figsize=(13,7));from scipy.special import expit
for ax,(s,(raw,a)) in zip(axs.flat,plan['sets'].items()):
 lo=min(d[raw].min(),a[0]);hi=max(a[2]*1.25,a[2]+(a[2]-a[0])*.1);x=np.linspace(lo,hi,400);l,c,u=a;f=expit(np.log(19)*(x-c)/np.where(x<c,c-l,u-c));ax.plot(x,f,color='#1b5674');ax.scatter(a,[.05,.5,.95],color='#bd6232',s=25);ax.axhline(.5,color='grey',lw=.7,ls='--');ax.set(title=s,xlabel={'FD':'Bank credit / GDP (%)','RQ':'WGI estimate units','GDP':'GDP per capita (constant 2015 USD)','Ypc':'Annual commitments / resident\n(constant 2022 USD)','Yabs':'Annual commitments\n(million constant 2022 USD)'}[s],ylabel='Membership',ylim=(0,1));ax.grid(alpha=.15)
axs.flat[-1].axis('off');fig.suptitle('Approved direct calibration — fixed raw-unit anchors');fig.tight_layout();save(fig,'01_calibration_curves')
fig,axs=plt.subplots(1,2,figsize=(12,4.6),sharey=True)
for ax,y in zip(axs,['Ypc','Yabs']):
 z=pt[pt.scenario=='primary_'+y];xx=np.arange(8);ax.bar(xx-.18,z.consistency,.36,label='Consistency',color='#28617a');ax.bar(xx+.18,z.PRI,.36,label='PRI',color='#c0804d');ax.axhline(.8,color='#28617a',ls='--',lw=.9);ax.axhline(.65,color='#c0804d',ls=':',lw=1);ax.set_xticks(xx,z.corner);ax.set(title=y+' — n: '+', '.join(map(str,z.n)),xlabel='FD / RQ / GDP corner',ylim=(0,1.04));ax.legend(loc='upper left',fontsize=8)
axs[0].set_ylabel('Set fit');fig.suptitle('Primary truth tables: thresholds applied jointly');fig.tight_layout();save(fig,'02_truth_table_fit')
fig,ax=plt.subplots(figsize=(7.5,4.5));z=pn[pn.scenario=='primary_Ypc'];ax.barh(z.condition,z.consistency,color='#28617a');ax.axvline(.9,color='#b65e3f',ls='--',label='Necessity screen .90');ax.set(xlim=(0,1),xlabel='Necessity consistency (screen = .90 dashed line)',title='No condition or absence passes primary necessity screen');save(fig,'03_primary_necessity')
fig,axs=plt.subplots(1,2,figsize=(11,4.5))
for ax,p in zip(axs,['RQ','FD']):
 for g,color in zip(groups,['#8b5b92','#4e83a7','#72a089','#b07848']):
  z=d[d.income_group_2018==g];ax.scatter(z[p],z.Yabs,s=26,label=g,color=color,alpha=.85)
 ax.plot([0,1],[0,1],ls='--',c='grey',lw=.8);ax.axhline(.5,lw=.5,c='grey');ax.axvline(.5,lw=.5,c='grey');ax.set(xlim=(0,1.02),ylim=(0,1.02),xlabel=p+' membership',ylabel='Yabs membership',title='Alternative volume pathway: '+p)
axs[1].legend(fontsize=8,loc='lower right');fig.tight_layout();save(fig,'04_alternative_pathway_case_plots')
fig,axs=plt.subplots(1,2,figsize=(11,4.4));cov=coverage[coverage.income_group_2018.isin(groups)].set_index('income_group_2018').loc[groups];axs[0].bar(range(4),cov.case_retention_share*100,color='#28617a');axs[0].set_xticks(range(4),['Low','Lower-middle','Upper-middle','High']);axs[0].set(ylabel='Eligible share of economy frame (%)',ylim=(0,100),title='Source-selection limitation');
for i,r in enumerate(cov.itertuples()):axs[0].text(i,r.case_retention_share*100+2,f'{r.eligible_cases}/{r.economies}',ha='center',fontsize=9)
matrix=cg.pivot(index='income_group',columns='corner',values='cases').loc[groups,corner];im=axs[1].imshow(matrix,aspect='auto',cmap='Blues');axs[1].set_xticks(range(8),corner,rotation=45);axs[1].set_yticks(range(4),['Low','Lower-middle','Upper-middle','High']);axs[1].set(title='Countries by condition corner')
for i in range(4):
 for j in range(8):axs[1].text(j,i,str(matrix.iloc[i,j]),ha='center',va='center',color='white' if matrix.iloc[i,j]>30 else '#203047')
fig.tight_layout();save(fig,'05_income_group_support')
fig,ax=plt.subplots(figsize=(9,4.5));za=aug[(aug.outcome=='Yabs')&(aug.kind=='leave_one_out')]
for same,color,label in [(True,'#28617a','Same complex expression (105 deletions)'),(False,'#bd6232','Changed complex expression (6 deletions)')]:
 zz=za[za.same_expression_as_primary==same];ax.scatter(zz.solution_coverage,zz.solution_consistency,color=color,s=25,alpha=.6,label=label)
pr=ps[(ps.scenario=='primary_Yabs')&(ps.type=='complex')].iloc[0];ax.scatter([pr.coverage],[pr.consistency],marker='*',s=180,color='black',label='Primary absolute-volume solution');ax.set(xlabel='Solution coverage',ylabel='Solution consistency',title='Country-deletion checks — absolute-volume outcome');ax.legend(fontsize=8);save(fig,'06_country_deletion_sensitivity')
# Inspectable source-selection contrasts, not outcome-based regression or formal selection correction.
sel=[]
for g in groups:
 for eligible in [True,False]:
  z=allcases[(allcases.income_group_2018==g)&(allcases.eligible==eligible)]
  for col in ['bank_depth_mean_2018_2020','RQ_mean_2018_2020','real_GDP_pc_mean_2018_2020']:
   sel.append(dict(income_group=g,eligible=eligible,variable=col,economies=len(z),observed=int(z[col].notna().sum()),median=z[col].median()))
selection=pd.DataFrame(sel);selection.to_csv(EV/'selection_observed_comparison.csv',index=False)
lo=aug[aug.kind=='leave_one_out'];freqs=aug.groupby(['outcome','kind','solution_expression'],dropna=False).size().reset_index(name='scenarios');freqs.to_csv(EV/'solution_frequency_by_sensitivity.csv',index=False)
def table(df,cols=None):return (df[cols] if cols else df).to_markdown(index=False,floatfmt='.4f')
hi=d[d.income_group_2018=='High income'];highnames=', '.join(hi.country_name)
primaryalt=ps[(ps.scenario=='primary_Yabs')&(ps.type=='complex')].iloc[0]
summary={'scenarios':len(meta),'country_deletions':int((meta.kind=='leave_one_out').sum()),'R_Python_checks_passed':bool(validation.passed.all()),'primary_Ypc_solutions':0,'alternative_expression':primaryalt.expression,'alternative_consistency':primaryalt.consistency,'alternative_PRI':primaryalt.PRI,'alternative_coverage':primaryalt.coverage,'abs_loo_same_expression':int(za.same_expression_as_primary.sum()),'abs_loo_count':len(za),'abs_loo_consistency_range':[float(za.solution_consistency.min()),float(za.solution_consistency.max())],'abs_loo_coverage_range':[float(za.solution_coverage.min()),float(za.solution_coverage.max())]}
(EV/'results_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
report=f'''# Approved fsQCA results — renewable energy public-financing commitments

## Decision and primary result

The approved primary analysis yields **no sufficient configuration for high per-capita selected-public-financing commitments**, and no individual condition or its absence passes the 0.90 necessity screen. This null configurational result is substantive evidence about the limits of this three-condition design, not a reason to retune its anchors. The alternative absolute-volume outcome yields **FD + RQ** in both complex and parsimonious minimization (Boolean OR, not an interaction), with consistency {primaryalt.consistency:.4f}, PRI {primaryalt.PRI:.4f}, coverage {primaryalt.coverage:.4f}. It describes a different financing construct and does not replace the primary result.

## Cases, timing, provenance and calibration

There are 111 unique country/economy cases: 24 low, 41 lower-middle, 35 upper-middle and 11 high income under historical 2018 classifications. All conditions require complete observed 2018–2020 values; both outcomes require complete recorded 2021–2023 commitments. Arithmetic annual means are used; per-capita outcome is the mean of annual commitments divided by that year's population, not cumulative commitments or a ratio of sums. Data include 87 cases whose historical income category stays constant through 2023. No cases, missing observations or vintages were appended; no missing financing was made zero. Published numerical zeros are retained with their original rounding flags. WDI reporting flags do not certify wholly unestimated publisher inputs. No IMF or FinData index was substituted for bank-depth data.

Financing comes exclusively from the frozen, independently PDF-validated IRENA Renewable Energy Statistics 2025: selected public renewable-energy commitments reported through the selected development-finance institutions and OECD DAC channel. It is not total national energy investment, disbursements, private finance, fiscal spending or verified absence of all other financing. The archived API remains unavailable; existing PDF reconciliation provides provisional source validity. The institutional perimeter, incomplete recipient reporting and unresolved 2025/2026 revisions (including UAE 2022) remain material limits. No 2026 financing values enter this analysis. Bank private credit/GDP measures intermediation depth, not efficiency or overall financial development. WGI estimates contain measurement uncertainty; three-year averaging does not eliminate it or establish regulatory adequacy at zero.

Approved anchors are researcher-specified materiality/reference judgments, not official policy thresholds or empirically estimated causal breakpoints:

{table(anchors,['set','full_out','crossover','full_in','members_above_crossover','exact_crossovers'])}

For l<c<u, membership is logistic(log(19)*(x-c)/(c-l)) below c and logistic(log(19)*(x-c)/(u-c)) at/above c. Anchors map to .05/.50/.95. No epsilon or tail clipping is applied. Custom R `plogis` and independent SciPy `expit` agree within {max(freeze['Python_R_max_errors'].values()):.3g}. QCA::calibrate differs in extreme GDP and absolute-outcome tails because it clips memberships <.0001 or >.9999 to 0/1; maximum differences {discrep.QCA_max_difference.max():.7f}. We retain the explicitly approved exact function and disclose, rather than silently adopt, that package convention. No primary score equals .5. Source and calibrated SHA256 are `{freeze['source_sha256']}` and `{freeze['calibrated_sha256']}`. The calibrated hash was frozen at {freeze['frozen_utc']}, before any truth-table calculations.

## Necessity

{table(pn,['scenario','condition','consistency','coverage','RoN','screen_pass'])}

Necessity consistency is sum min(X,Y)/sum Y; necessity coverage sum min(X,Y)/sum X. RoN=(N−sum X)/(N−sum min(X,Y)) guards against trivial supersets but is not a causal identification statistic. None passes .90 for either outcome. Contradiction mass and named high-outcome cases lacking each candidate condition are preserved in the necessity workbook. No necessity claim is supported; changing the local screen is reported only as sensitivity.

## Complete truth tables and sufficiency

{table(pt,['scenario','corner','n','consistency','PRI','OUT','contradictory_n'])}

Corner digits are FD/RQ/GDP, 1 for presence and 0 for absence. Configuration membership is the minimum of its literals; frequency counts only memberships >.5. Fuzzy fit statistics use all retained countries, not just corner owners. Exact .5 cases, if encountered in sensitivity, receive no corner and remain in fit sums. OUT=1 requires frequency≥3, consistency≥.80 and PRI≥.65 jointly. OUT=? denotes insufficient frequency, not an assumed positive or negative outcome. All eight primary corners are observed; 001 has only two cases and is a frequency remainder. Thus scarcity is mainly highly uneven support, not wholly unobserved corners. Primary per-capita PRI values are all below .65, including the corner with consistency >.80. High consistency alone would conceal contradictory outcome overlap.

Consistency=sum min(X,Y)/sum X; coverage=sum min(X,Y)/sum Y; PRI=[sum min(X,Y)−sum min(X,Y,1−Y)]/[sum X−sum min(X,Y,1−Y)]. PRI is not computed by subtracting sum min(X,1−Y). Country owners and contradictory cases (corner>.5 and Y<.5) are listed for every row.

R QCA::truthTable was called with complete=TRUE, n.cut=3, incl.cut=c(.80,.80), pri.cut=.65: matching inclusion cutpoints prevent an additional ambiguous C state. QCA::minimize uses include='' for complex and include='?' for parsimonious, row.dom=FALSE and all.sol=TRUE. The independent Python implementation enumerates all 27 possible three-condition cubes, prime implicants and every minimum-term cover. Complex minimization excludes frequency remainders; parsimonious permits them but excludes supported negative rows. All competing minimal covers are retained. No minimization was forced when positive configurations were absent. Intermediate solutions are withheld because credible opposing theoretical directions for selected public-recipient finance preclude defensible directional counterfactuals. Therefore no core/peripheral labels are claimed merely from complex-versus-parsimonious agreement.

The absolute-volume complex and parsimonious solutions coincide and require no frequency-remainder counterfactual to cover their positive rows. Pathway membership uses fuzzy minimum and solution membership fuzzy maximum. Raw coverage is pathway–outcome overlap; unique coverage is solution–outcome overlap minus overlap explained by the union of the other terms, divided by total outcome membership:

{table(pp[['type','pathway','consistency','PRI','raw_coverage','unique_coverage']])}

`FD + RQ` means either high bank depth or relatively high RQ in this selected set; it does not establish that depth and governance substitute causally, nor an FD×RQ interaction. Absolute amounts are sensitive to national scale, large projects, commitment timing and shared recipients/providers. All supporting and contradictory countries appear in the pathway workbook, including 3 RQ owners and 7 FD owners below outcome crossover. These are fallible descriptive pathways, not universal rules.

## Prespecified robustness

All {len(meta)} scenarios are preserved, including {len(lo)} leave-one-country-out runs (111 per outcome). Checks are one-factor-at-a-time, not a threshold search. Endpoints remain fixed when approved crossovers change. Frequency 2/3/5, local consistency .75/.80/.85, PRI .60/.65/.70 and necessity .85/.90/.95 are assessed. The user authorized criterion sensitivity; those local numeric alternatives were operationalized and frozen before calibration, not retrospectively described as old numeric preregistrations. Old approved protocol also permits absolute crossover 10/50, log-distance depth/GDP, log1p currency-unit-distance outcomes and the per-capita ratio-of-sums alternative; all are separately labelled. Fixed near-crossover bands are dropped individually and jointly; these are sensitivity samples, never replacements for the approved 111 cases.

No primary per-capita configuration appears in any permitted sensitivity scenario, including all country deletions. Absolute-volume solution expressions vary across calibrations and thresholds; complete scenario-level positive rows, competing covers, fit, failures and actual country support are supplied. In country deletions, {summary['abs_loo_same_expression']}/{summary['abs_loo_count']} retain the exact complex expression FD + RQ; consistency ranges {za.solution_consistency.min():.4f}–{za.solution_consistency.max():.4f}, coverage {za.solution_coverage.min():.4f}–{za.solution_coverage.max():.4f}. Country deletion changes frequencies and thus can alter Boolean expressions even where fit changes little. Robustness is described across the full menu; no favourable variant is promoted to primary. No outcome-negation minimization or unapproved extra-condition search was run.

## Income context and selection

{table(grp,['income_group','cases','high_Ypc','high_Yabs','observed_configurations','stable_income_2018_2023'])}

The 11 high-income cases are {highnames}. They are 11/80 (13.75%) of the reference-frame high-income economies; 68/80 lack a complete financing window, overlapping other exclusions. The reference frame includes territories/SARs and must not be called 80 sovereign states. This is severe source/recipient selection, not a random sample. The selected high-income cases differ from missing economies, whose financing cannot be assigned a zero. Descriptive available-condition medians and observed counts for retained/excluded economies are in audit_evidence/selection_observed_comparison.csv, without formal correction or missing-not-at-random assumptions. The subgroup tables retain all contexts but no separate confirmatory four-group minimizations were run. Configuration counts of a few high-income cases do not establish high-income mechanisms or support extrapolation to all high-income financing systems.

## Independent verification and remaining limits

R 4.6.1, QCA 3.25, admisc .41, declared .27 and venn 1.13 were used. All {len(validation)} numerical/discrete R–Python validation checks pass (numerical tolerance 1e-10; discrete exact). Membership verification, all {len(tt)} truth-table rows, {len(nec)} necessity rows, all nonempty solutions and pathway fit/unique coverage are independently reconciled. Package archives, URLs, hashes, stdout/stderr and session information are retained. CRAN index did not supply MD5 fields for these downloads; SHA256 provenance was recorded, not falsely described as externally MD5-verified. R emitted locale-startup warnings; package calculations completed and numerical checks passed. Author-source downloads are methodological references, not a claim their branch equals the installed release.

Calibration uncertainty is assessed through declared sensitivity only, not statistical confidence intervals. Cases are cross-nationally dependent through providers, projects and regional shocks; one country per case prevents repeated-year pseudoreplication but does not prove stochastic independence. Temporally ordered windows reduce contemporaneous overlap, not reverse causality, selection, omitted factors or financing anticipation. COVID-era conditions and three-year project timing limit temporal generalization. There is no causal identification, universal necessity/sufficiency, significance test, or solution validation from agreement with panel signs. Preserve the original inconclusive panel findings and present this fsQCA as a transparent complementary, exploratory study.

## Files and reproducibility

All requested workbooks, calibrated CSV, six PNG/PDF figures and scripts are in this results directory. `scripts/README.md` states execution order, dependencies, freeze gates and equations. `audit_evidence/independent_validation.csv`, `calibration_frozen_receipt.json` and `final_validation.json` provide numerical, temporal and file-preservation receipts. Prior panel results and preparation/source files were not changed.
'''
report+='''\n## Additional identity, sensitivity and display audit\n\nThe 111 cases comprise 109 UN member states, one UN observer State and Kosovo as a separate World Bank economy (XKX retained as the project identifier, not represented as an assigned ISO3166 code). No aggregate regions are treated as cases. Country-deletion expression changes occur when Brazil, China or Türkiye is omitted (FD*~GDP + RQ), or Kazakhstan, Mexico or Romania is omitted (FD + RQ*~GDP). These deletions reduce three-case configurations below the fixed frequency threshold; they do not motivate case exclusion. Exact expressions are more fragile to RQ/GDP anchors and frequency than aggregate fit values.\n\nThe public R truthTable display masks fit for some zero-frequency sensitivity corners with a dash. Its original output is preserved. Full-precision fits for all corners are separately exported using the installed package C_truthTable routine and reconciled independently; missing display values were not assigned zero or silently treated as calculations.\n'''
(OUT/'fsqca_results_report.md').write_text(report,encoding='utf-8')
integrated=f'''# Panel and configurational evidence: integrated interpretation

The analyses address related but different questions. The panel models examine conditional within-country associations in 2014–2023 with country/year effects; fsQCA examines set relations among three-year baseline conditions (2018–2020) and subsequent financing (2021–2023) for 111 selected countries. A condition's fuzzy membership is not a panel coefficient or marginal effect.

The replicated panel evidence remains unchanged: financial depth is inconclusive; RQ is positive in OLS but sensitive in PPML; FD×RQ and income-group contrasts are inconclusive. RISE, FDI and access-deficit mechanism tests are inconclusive. The original benchmark remains 1,104 observations/137 clusters. No new regressions were run and no panel specification or reported coefficient was replaced.

The primary per-capita fsQCA finds no condition or absence meeting necessity .90 and no corner meeting frequency3, consistency.80, PRI.65 jointly. This persists through all prespecified checks. It supplies no evidence for the approved three-condition sufficient pathways to high recorded public-finance intensity. It neither proves no relationships exist nor validates an omitted mechanism.

The distinct absolute-volume outcome yields FD + RQ, consistency {primaryalt.consistency:.4f}, PRI {primaryalt.PRI:.4f}, coverage {primaryalt.coverage:.4f}, in both complex and parsimonious solutions. These are alternative routes in Boolean notation (OR), not a depth–governance interaction or evidence of causal substitutability. The full sensitivity record qualifies robustness; expression changes are disclosed, not screened away. Population scale and large transactions can elevate volume while leaving per-capita intensity low. Therefore the absolute result cannot be relabelled as a successful primary intensity model or confirmation of inconclusive panel interactions.

What fsQCA adds is an explicit audit of configurational sufficiency/necessity, case ownership, contradictions, outcome-definition dependence and counterfactual reliance. What it does not add is identification of causal effects, rehabilitation of insignificant hypotheses, or external validity to missing high-income economies. Only 11 of 80 high-income reference-frame economies enter; missing financing is not zero. Historical income groups remain context, not four independent confirmatory models.

The common outcome is selected IRENA2025 public commitments, not all renewable investment. Source selection, provider limits, archived-API unavailability, unresolved release revisions, WGI uncertainty and project timing remain shared limitations. Temporal ordering does not eliminate anticipation or reverse causality. No intermediate directional solution or core/peripheral causal label is justified.

For the article, retain panel estimates and prespecified inference as the primary associational evidence. Report fsQCA separately as exploratory configurational evidence, leading with its primary null sufficiency result and necessity failure. Present absolute volume only as the predefined alternative construct. Include the entire sensitivity menu and high-income coverage table. The mixed evidence motivates restraint and better measurements, not specification or calibration changes aimed at significance.
'''
(OUT/'panel_fsqca_integrated_interpretation.md').write_text(integrated,encoding='utf-8')
versions={'Python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__,'scipy':__import__('scipy').__version__,'matplotlib':matplotlib.__version__,'R':'4.6.1','QCA':'3.25'}
(SCRIPTS/'dependency_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
for name in ['qca_setup.py','qca_calibrate.R','qca_run_calibration.py','qca_analysis.R','qca_compute.py','qca_finish.py','qca_packages.py','qca_precision.py','qca_precision.R']:
 shutil.copy2(stage/name,SCRIPTS/name)
print(json.dumps(summary,indent=2));print(validation.to_string(index=False));print(freqs.to_string(index=False))
