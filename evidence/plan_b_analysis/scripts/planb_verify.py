from planb_setup import *
import numpy as np,pandas as pd,importlib.util
from scipy.linalg import qr,eigh
from scipy.optimize import minimize_scalar
from scipy import stats
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
spec=importlib.util.spec_from_file_location('frozen_engine',ROOT/'econometric_results/econometric_engine.py');eng=importlib.util.module_from_spec(spec);spec.loader.exec_module(eng)
plan=pd.read_csv(EV/'model_plan.csv');rco=pd.read_csv(EV/'R_coefficients.csv');rmeta=pd.read_csv(EV/'R_model_metadata.csv');rjoint=pd.read_csv(EV/'R_joint_tests.csv');rm=rmeta.set_index('model');baseX=['FD','RQ','lnGDP','lnPop','Trade','Inflation']
checks=[];coeffs=[];met=[];covs=[];joint=[];contrasts=[];allfits={};variancechecks=[]
class Independent(eng.Fit):
 def __init__(self,z,D,ratio=0):
  self.z=z.reset_index(drop=True);z=self.z;self.groups=[np.flatnonzero(z.country_code==g) for g in sorted(z.country_code.unique())];self.codes=sorted(z.country_code.unique());self.G=len(self.codes);self.N=len(z)
  names=list(D.columns);X=D.to_numpy().copy();y=z.y.to_numpy().copy()
  for ix in self.groups:
   c=(1/np.sqrt(1+ratio*len(ix))-1)/len(ix);X[ix]+=c*X[ix].sum(axis=0);y[ix]+=c*y[ix].sum()
  norm=np.linalg.norm(X,axis=0);Z=X/norm;self.Q,R=qr(Z,mode='economic');assert np.min(abs(np.diag(R)))>1e-10
  self.rank=Z.shape[1];self.beta_full=np.linalg.lstsq(Z,y,rcond=None)[0]/norm;self.e=y-X@self.beta_full;self.X=X;self.y=y
  focus=baseX if 'W_FD' not in names else ['W_'+c for c in baseX]+['B_'+c for c in baseX];idx=[names.index(c) for c in focus]
  self.names=focus;self.p=len(idx);self.beta=self.beta_full[idx];self.B=np.linalg.solve(Z.T@Z,Z.T)[idx].T/norm[idx]
  self.infer();self.condition=float(np.linalg.cond(Z));self.VIF=np.diag(np.linalg.inv(Z.T@Z))[idx]
def reml_profile(z,D):
 X=D.to_numpy();X=X/np.linalg.norm(X,axis=0);y=z.y.to_numpy();groups=[np.flatnonzero(z.country_code==g) for g in sorted(z.country_code.unique())];N,p=X.shape
 def objective(logratio):
  ratio=np.exp(logratio);XX=X.copy();yy=y.copy();ld=0
  for ix in groups:
   n=len(ix);c=(1/np.sqrt(1+ratio*n)-1)/n;XX[ix]+=c*X[ix].sum(axis=0);yy[ix]+=c*y[ix].sum();ld+=np.log1p(ratio*n)
  beta=np.linalg.lstsq(XX,yy,rcond=None)[0];e=yy-XX@beta;sign,ldx=np.linalg.slogdet(XX.T@XX);return ld+ldx+(N-p)*np.log(e@e/(N-p))
 opt=minimize_scalar(objective,bounds=(-18,18),method='bounded',options={'xatol':1e-8});return float(np.exp(opt.x)),float(opt.fun)
with threadpool_limits(limits=1):
 for m in plan.to_dict('records'):
  name=m['model'];z=pd.read_csv(EV/(name+'_data.csv'));D=pd.read_csv(EV/(name+'_design.csv'));ratio=rm.loc[name,'tau2']/rm.loc[name,'sigma2'] if m['method']=='CRE' else 0
  f=Independent(z,D,ratio);allfits[name]=f
  for r in f.coefficients():coeffs.append(dict(model=name,**r))
  met.append(dict(model=name,n=f.N,clusters=f.G,rank=f.rank,residual_df=f.N-f.rank,condition_index=f.condition,max_VIF=f.VIF.max(),working_cov_identity_error=f.working_model_covariance_error))
  for i,t in enumerate(f.names):
   for j,u in enumerate(f.names):covs.append(dict(model=name,term1=t,term2=u,covariance_CR2=f.cov[i,j]))
  if m['method']=='CRE':
   prof,obj=reml_profile(z,D);variancechecks.append(dict(model=name,R_ratio=ratio,Python_profile_REML_ratio=prof,relative_difference=abs(prof-ratio)/max(1,ratio)))
  if name=='annual_CRE_absolute':
   for vars in [baseX,['FD'],['RQ']]:joint.append(dict(model=name,test='within_equals_between_'+';'.join(vars),**f.joint([f.unit('W_'+v)-f.unit('B_'+v) for v in vars],'Within equals between')))
   for v in ['FD','RQ']:contrasts.append(dict(model=name,**f.contrast(f.unit('W_'+v)-f.unit('B_'+v),'W_'+v+' minus B_'+v)))
  print('PYTHON_VERIFIED',name,flush=True)
ours=pd.DataFrame(coeffs);comparison=ours.merge(rco,on=['model','term'],suffixes=('_Python','_R'));assert len(comparison)==len(ours)==len(rco)
for col in ['estimate','SE_CR2','Satterthwaite_df','p_value','CI95_low','CI95_high']:
 a=comparison[col+'_Python'];b=comparison[col+'_R'];error=float(np.max(abs(a-b)));scaled=float(np.max(abs(a-b)/(1+abs(b))));checks.append(dict(check='R–Python '+col,max_absolute_difference=error,max_scaled_difference=scaled,tolerance=2e-6,passed=scaled<2e-6))
pj=pd.DataFrame(joint);jc=pj.merge(rjoint,on=['model','test'],suffixes=('_Python','_R'))
for a,b in [('HTZ_F','HTZ_F'),('HTZ_df2','HTZ_df2'),('HTZ_p','p_value')]:
 error=float(np.max(abs(jc[a+'_Python' if a in rjoint.columns else a]-jc[b+'_R' if b in pj.columns else b])));checks.append(dict(check='R–Python joint '+a,max_absolute_difference=error,tolerance=2e-5,passed=error<2e-5))
for name in ['annual','balanced','lagged']:
 f=allfits[name+'_TWFE'+('_absolute' if name=='annual' else '')];c=allfits[name+'_CRE'+('_absolute' if name=='annual' else '')];err=max(abs(f.beta[f.names.index(x)]-c.beta[c.names.index('W_'+x)]) for x in baseX);checks.append(dict(check=name+' matched TWFE–CRE within coefficient identity',max_absolute_difference=err,tolerance=1e-7,passed=err<1e-7))
old=pd.read_csv(ROOT/'econometric_results/audit_evidence/coefficients.csv');old=old[old.model=='Primary Baseline'];mapping=dict(FD='FD',RQ='Governance',lnGDP='log_real_GDP_pc',lnPop='log_population',Trade='trade_GDP',Inflation='inflation')
for x,legacy in mapping.items():
 new=ours[(ours.model=='annual_TWFE_absolute')&(ours.term==x)].iloc[0];oldrow=old[old.term==legacy].iloc[0]
 err=max(abs(new[k]-oldrow[k]) for k in ['estimate','SE_CR2','Satterthwaite_df','p_value','CI95_low','CI95_high']);checks.append(dict(check='Frozen original benchmark '+x,max_absolute_difference=err,tolerance=1e-6,passed=err<1e-6))
vchecks=pd.DataFrame(variancechecks);vchecks['passed']=vchecks.relative_difference<2e-4
checks.append(dict(check='Independent profile REML variance ratios',max_absolute_difference=float(vchecks.relative_difference.max()),tolerance=2e-4,passed=bool(vchecks.passed.all())))
# Apply the frozen eight-test Holm family, without selecting a specification.
ours['p_Holm_primary8']=np.nan;family=((ours.model=='annual_CRE_absolute')&ours.term.isin(['W_FD','B_FD','W_RQ','B_RQ'])) | (ours.model.isin(['window2_TWFE_absolute','window2_TWFE_per_capita'])&ours.term.isin(['FD','RQ']));assert family.sum()==8;ours.loc[family,'p_Holm_primary8']=multipletests(ours.loc[family,'p_value'],method='holm')[1]
ct=pd.DataFrame(contrasts);ct['p_Holm_two_differences']=multipletests(ct.p_value,method='holm')[1]
pd.DataFrame(checks).to_csv(EV/'independent_validation.csv',index=False);comparison.to_csv(EV/'independent_coefficient_comparison.csv',index=False);jc.to_csv(EV/'independent_joint_comparison.csv',index=False);vchecks.to_csv(EV/'independent_REML_variance_comparison.csv',index=False)
assert pd.DataFrame(checks).passed.all(),pd.DataFrame(checks).to_string(index=False)
for name,df in [('coefficients',ours),('Python_model_metadata',pd.DataFrame(met)),('Python_CR2_covariances',pd.DataFrame(covs)),('within_between_contrasts',ct)]:df.to_csv(EV/(name+'.csv'),index=False)
# Country deletions are diagnostics, not new preferred fits. CRE variance ratio fixed.
loco=[]
with threadpool_limits(limits=1):
 for name in ['annual_CRE_absolute','window2_TWFE_absolute']:
  z=pd.read_csv(EV/(name+'_data.csv'));D=pd.read_csv(EV/(name+'_design.csv'));f=allfits[name];ratio=rm.loc[name,'tau2']/rm.loc[name,'sigma2'] if name.startswith('annual') else 0
  focus=['W_FD','B_FD','W_RQ','B_RQ'] if name.startswith('annual') else ['FD','RQ']
  for code in sorted(z.country_code.unique()):
   keep=z.country_code!=code;zz=z[keep].reset_index(drop=True);XX=D[keep].reset_index(drop=True);names=list(XX.columns);A=XX.to_numpy();y=zz.y.to_numpy().copy()
   for ix in [np.flatnonzero(zz.country_code==g) for g in sorted(zz.country_code.unique())]:
    c=(1/np.sqrt(1+ratio*len(ix))-1)/len(ix);A[ix]+=c*A[ix].sum(axis=0);y[ix]+=c*y[ix].sum()
   norm=np.linalg.norm(A,axis=0);active=norm>1e-10;beta=np.linalg.lstsq(A[:,active]/norm[active],y,rcond=None)[0]/norm[active];bn=dict(zip(np.array(names)[active],beta))
   for x in focus:loco.append(dict(model=name,omitted_country=code,term=x,estimate=bn[x],full_estimate=f.beta[f.names.index(x)],change=bn[x]-f.beta[f.names.index(x)],variance_treatment='Fixed full-sample REML variance ratio; conditional FGLS diagnostic' if ratio else 'OLS refitted; country/window effects retained'))
pd.DataFrame(loco).to_csv(EV/'country_deletion_diagnostics.csv',index=False)
print('ALL INDEPENDENT CHECKS PASSED',flush=True);print(ours[family].to_string(index=False));print(ct.to_string(index=False))
