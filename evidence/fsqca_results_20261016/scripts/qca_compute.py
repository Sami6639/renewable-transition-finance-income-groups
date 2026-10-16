from qca_setup import *
import numpy as np,pandas as pd,itertools,subprocess,re,platform
from scipy.special import expit
stage=Path(__file__).parent
plan=json.loads((EV/'approved_execution_plan_frozen.json').read_text());receipt=json.loads((EV/'calibration_frozen_receipt.json').read_text());assert sha(OUT/'fsqca_calibrated_dataset.csv')==receipt['calibrated_sha256']
d=pd.read_csv(OUT/'fsqca_calibrated_dataset.csv');sets=plan['sets'];conditions=['FD','RQ','GDP'];corners=list(itertools.product([0,1],repeat=3))
def calibrate(x,a,shape='raw'):
 x=np.asarray(x);a=np.array(a,dtype=float)
 if shape=='log':x=np.log(x);a=np.log(a)
 elif shape=='log1p':x=np.log1p(x);a=np.log1p(a)
 l,c,u=a;return expit(np.log(19)*(x-c)/np.where(x<c,c-l,u-c))
scenarios=[];frames=[]
def add(name,outcome,kind='primary',change='',n=3,inc=.8,pri=.65,nec=.9,shift=None,shape=None,keep=None,denominator=False):
 z=d[['country_code','country_name','income_group_2018','FD','RQ','GDP']].copy();z['Y']=d[outcome];z['scenario']=name
 if shift:
  s,c=shift;a=list(sets[s][1]);a[1]=c;z[s if s in conditions else 'Y']=calibrate(d[sets[s][0]],a)
 if shape:
  s,sh=shape;z[s if s in conditions else 'Y']=calibrate(d[sets[s][0]],sets[s][1],sh)
 if denominator:z['Y']=calibrate(d.financing_ratio_of_sums_per_capita,sets['Ypc'][1])
 if keep is not None:z=z.loc[keep].copy()
 scenarios.append(dict(scenario=name,outcome=outcome,kind=kind,change=change,cases=len(z),frequency=n,consistency=inc,PRI=pri,necessity=nec));frames.append(z)
for y in ['Ypc','Yabs']:
 add('primary_'+y,y)
 for s,cs in {'FD':[40,60],'RQ':[-.25,.25],'GDP':[6000,10000],y:([5,20] if y=='Ypc' else [10,50])}.items():
  for c in cs:add(f'{y}_cross_{s}_{c}',y,'crossover',f'{s} crossover {c}; fixed endpoints',shift=(s,c))
 for n in [2,5]:add(f'{y}_frequency_{n}',y,'frequency',str(n),n=n)
 for ic in [.75,.85]:add(f'{y}_consistency_{ic}',y,'criterion',f'consistency {ic}',inc=ic)
 for p in [.60,.70]:add(f'{y}_PRI_{p}',y,'criterion',f'PRI {p}',pri=p)
 for nec in [.85,.95]:add(f'{y}_necessity_{nec}',y,'necessity',f'necessity {nec}',nec=nec)
 for s,sh in [('FD','log'),('GDP','log'),(y,'log1p')]:add(f'{y}_shape_{s}',y,'shape',s+' '+sh,shape=(s,sh))
 if y=='Ypc':add('Ypc_ratio_of_sums',y,'denominator','Ratio of sums rather than mean annual ratios',denominator=True)
 bands={'FD':5,'RQ':.1,'GDP':800,y:1 if y=='Ypc' else 2.5};near={s:(abs(d[sets[s][0]]-sets[s][1][1])<=b) for s,b in bands.items()}
 for s,mask in near.items():add(f'{y}_near_drop_{s}',y,'near_crossover',f'Drop raw-unit band {s} +/-{bands[s]}',keep=~mask)
 add(f'{y}_near_drop_conditions',y,'near_crossover','Drop union of condition bands',keep=~np.logical_or.reduce([near[s] for s in conditions]))
 add(f'{y}_near_drop_all',y,'near_crossover','Drop union of conditions and outcome bands',keep=~np.logical_or.reduce(list(near.values())))
 for code in d.country_code:add(f'{y}_leave_{code}',y,'leave_one_out',code,keep=d.country_code!=code)
meta=pd.DataFrame(scenarios);members=pd.concat(frames,ignore_index=True);meta.to_csv(EV/'scenario_metadata.csv',index=False);members.to_csv(EV/'scenario_memberships.csv',index=False)
print('Frozen scenarios',len(meta),'membership rows',len(members),flush=True)
def fit(x,y):
 a=np.minimum(x,y).sum();b=np.minimum(np.minimum(x,y),1-y).sum();return dict(consistency=a/x.sum(),PRI=(a-b)/(x.sum()-b),coverage=a/y.sum())
def cubemem(cube,z):return np.minimum.reduce([z[c].values if b else 1-z[c].values for c,b in zip(conditions,cube) if b!=-1]) if cube!=(-1,-1,-1) else np.ones(len(z))
def cubetext(cube):return '*'.join(('' if b else '~')+c for c,b in zip(conditions,cube) if b!=-1) or '1'
def cubecorners(cube):return {b for b in corners if all(x==-1 or x==v for x,v in zip(cube,b))}
def minimize(pos,allow):
 if not pos:return []
 cand=[c for c in itertools.product([-1,0,1],repeat=3) if cubecorners(c)<=allow and cubecorners(c)&pos]
 primes=[c for c in cand if not any(cubecorners(c)<cubecorners(o) for o in cand)]
 for k in range(1,len(primes)+1):
  covers=[cs for cs in itertools.combinations(primes,k) if set.union(*(cubecorners(c) for c in cs))>=pos]
  if covers:return covers
 return []
tts=[];necs=[];solutions=[];paths=[];assumptions=[];case_support=[]
for m,z in zip(scenarios,frames):
 name=m['scenario'];y=z.Y.values;pos=set();rem=set()
 for b in corners:
  x=cubemem(b,z);f=fit(x,y);n=int((x>.5).sum());state='?' if n<m['frequency'] else '1' if f['consistency']>=m['consistency'] and f['PRI']>=m['PRI'] else '0'
  if state=='1':pos.add(b)
  if state=='?':rem.add(b)
  owner=x>.5;codes=z.country_code[owner];contr=owner&(y<.5)
  tts.append(dict(scenario=name,corner=''.join(map(str,b)),n=n,OUT=state,consistency=f['consistency'],PRI=f['PRI'],coverage=f['coverage'],countries=';'.join(codes),contradictory_cases=';'.join(z.country_code[contr]),contradictory_n=int(contr.sum()),unobserved=n==0,insufficient_frequency=n<m['frequency']))
 for c in conditions:
  for neg in [False,True]:
   x=1-z[c].values if neg else z[c].values;a=np.minimum(x,y).sum();cons=a/y.sum();contr=(y>.5)&(x<.5)
   necs.append(dict(scenario=name,condition=('~' if neg else '')+c,consistency=cons,coverage=a/x.sum(),RoN=(1-x).sum()/(1-np.minimum(x,y)).sum(),screen_pass=cons>=m['necessity'],contradiction_mass=float(np.maximum(y-x,0).sum()),contradictory_cases=';'.join(z.country_code[contr]),supporting_cases=';'.join(z.country_code[(x>.5)&(y>.5)]),high_outcome_cases=int((y>.5).sum())))
 for typ,allow in [('complex',pos),('parsimonious',pos|rem)]:
  covers=minimize(pos,allow)
  if not covers:solutions.append(dict(scenario=name,type=typ,model=0,expression='',consistency=np.nan,PRI=np.nan,coverage=np.nan,status='No sufficient configuration'));continue
  for mi,cubes in enumerate(covers,1):
   xs=[cubemem(c,z) for c in cubes];union=np.maximum.reduce(xs);f=fit(union,y)
   expr=' + '.join(sorted(cubetext(c) for c in cubes));used=set.union(*(cubecorners(c) for c in cubes))&rem
   solutions.append(dict(scenario=name,type=typ,model=mi,expression=expr,**f,terms=len(cubes),literals=sum(sum(b!=-1 for b in c) for c in cubes),counterfactual_corners=';'.join(''.join(map(str,b)) for b in sorted(used)),union_pass=f['consistency']>=m['consistency'] and f['PRI']>=m['PRI'],status='Boolean solution, not causal identification'))
   for pi,(cube,x) in enumerate(zip(cubes,xs)):
    others=np.maximum.reduce([v for k,v in enumerate(xs) if k!=pi]) if len(xs)>1 else np.zeros(len(z));unique=(np.minimum(union,y).sum()-np.minimum(others,y).sum())/y.sum();pf=fit(x,y)
    paths.append(dict(scenario=name,type=typ,model=mi,pathway=cubetext(cube),consistency=pf['consistency'],PRI=pf['PRI'],raw_coverage=pf['coverage'],unique_coverage=unique,countries=';'.join(z.country_code[x>.5]),supporting_high_outcome=';'.join(z.country_code[(x>.5)&(y>.5)]),contradictions=';'.join(z.country_code[(x>.5)&(y<.5)]),core_peripheral='Not identified: no justified intermediate solution'))
    if name.startswith('primary'):
     for idx,case in z.iterrows():
      loc=z.index.get_loc(idx);case_support.append(dict(scenario=name,type=typ,model=mi,pathway=cubetext(cube),country_code=case.country_code,country_name=case.country_name,income_group=case.income_group_2018,pathway_membership=x[loc],outcome_membership=y[loc],pathway_owner=bool(x[loc]>.5),supports_high_outcome=bool(x[loc]>.5 and y[loc]>.5),contradiction=bool(x[loc]>.5 and y[loc]<.5)))
for name,rows in [('Python_truth_tables',tts),('Python_necessity',necs),('Python_solutions',solutions),('Python_pathways',paths),('primary_pathway_case_support',case_support)]:pd.DataFrame(rows).to_csv(EV/(name+'.csv'),index=False)
shutil.copy2(stage/'qca_analysis.R',SCRIPTS/'qca_analysis.R')
r=Path(__import__('os').environ.get('RSCRIPT','Rscript'))
p=subprocess.run([str(r),'--vanilla',str(SCRIPTS/'qca_analysis.R'),str(ROOT)],capture_output=True,text=True,encoding='utf-8',errors='replace')
(EV/'R_analysis_stdout.txt').write_text(p.stdout,encoding='utf-8');(EV/'R_analysis_stderr.txt').write_text(p.stderr,encoding='utf-8');print('R return code',p.returncode,flush=True);print(p.stderr[-3000:]);assert p.returncode==0
for name,body in re.findall(r'BEGIN_EXPORT:([^\r\n]+)\r?\n(.*?)\r?\nEND_EXPORT',p.stdout,re.S):(EV/name).write_text(body+'\n',encoding='utf-8')
print('R and independent Python analysis complete',flush=True)
