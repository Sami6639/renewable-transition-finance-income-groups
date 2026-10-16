"""OLS, generalized CR2/identity-target inference, HTZ and restricted wild bootstrap.
No fabricated observations. Singular cluster directions from nested country FE are
handled by the Moore-Penrose inverse on the residual subspace.
"""
import numpy as np,pandas as pd
from scipy.linalg import qr,eigh
from scipy import stats
from statsmodels.stats.multitest import multipletests
G=['Low income','Lower middle income','Upper middle income','High income']
Y='irena2025_pdf';FD='FD.AST.PRVT.GD.ZS';RQ='regulatory_quality_estimate';MAC=['gdp_per_capita_constant_2015_usd','population','trade_pct_gdp','inflation_cpi_pct']
def core(z,ppml=False):
 while len(z):
  n=len(z);ok=z.groupby('country_code').year.transform('size').ge(2)&z.groupby('year').country_code.transform('size').ge(2)
  if ppml:ok &=z.groupby('country_code')[Y].transform('sum').gt(0)&z.groupby('year')[Y].transform('sum').gt(0)
  z=z[ok].copy()
  if len(z)==n:break
 return z.reset_index(drop=True)
def eligible(d,fd=FD,gov=RQ,pre=''):
 cols=[pre+fd,pre+gov]+[pre+c for c in MAC]
 return d[np.isfinite(d[cols+[Y]].astype(float)).all(axis=1)&d.group.isin(G)&d[pre+MAC[0]].gt(0)&d[pre+'population'].gt(0)&d[Y].ge(0)].copy()
def basis(X):
 Q,R,piv=qr(X,mode='economic',pivoting=True);rank=int((np.abs(np.diag(R))>max(X.shape)*np.finfo(float).eps*abs(R[0,0])).sum());return Q[:,:rank],rank,piv[:rank]
def variables(z,fd,gov,pre,block,centers=None):
 means=centers or (z[pre+fd].mean(),z[pre+gov].mean());v=pd.DataFrame({'FD':z[pre+fd]-means[0],'Governance':z[pre+gov]-means[1],'log_real_GDP_pc':np.log(z[pre+MAC[0]]),'log_population':np.log(z[pre+'population']),'trade_GDP':z[pre+'trade_pct_gdp'],'inflation':z[pre+'inflation_cpi_pct']})
 if block in ['RQ interaction','Combined']:v['FD_x_Governance']=v.FD*v.Governance
 if block in ['Income interaction','Combined']:
  for g in G[1:]:v['group_'+g]=z.group.eq(g).astype(float);v['FD_x_'+g]=v.FD*z.group.eq(g)
 return v,means
class Fit:
 def __init__(self,z,fd=FD,gov=RQ,pre='',block='Baseline',outcome='log1p',centers=None,inference=True,unit=1):
  self.z=z.reset_index(drop=True);z=self.z;self.fd=fd;self.gov=gov;self.pre=pre;self.block=block
  D=np.column_stack([pd.get_dummies(z.country_code,dtype=float),pd.get_dummies(z.year,drop_first=True,dtype=float)]);QD,rd,piv=basis(D);self.D=D[:,piv];self.QD=QD
  v,self.centers=variables(z,fd,gov,pre,block,centers);V=v.to_numpy();W=V-QD@(QD.T@V);W-=QD@(QD.T@W);norm=np.linalg.norm(W,axis=0);active=norm>1e-9*np.maximum(np.linalg.norm(V,axis=0),1);self.names=list(v.columns[active]);self.absorbed=list(v.columns[~active]);self.V=V[:,active];self.norm=norm[active];Z=W[:,active]/self.norm
  QS,rs,_=basis(Z);assert rs==Z.shape[1],('unidentified slope block',self.names);self.Z=Z;self.Q=np.column_stack([QD,QS]);self.rank=rd+rs;self.p=rs;self.N=len(z);self.G=z.country_code.nunique();self.year_count=z.year.nunique()
  self.B=(Z@np.linalg.inv(Z.T@Z))/self.norm;self.y=np.log1p(z[Y].to_numpy()*unit) if outcome=='log1p' else np.log(z[Y].to_numpy()) if outcome=='positive_log' else z[Y].to_numpy()
  self.beta=self.B.T@self.y;self.fitted=self.Q@(self.Q.T@self.y);self.e=self.y-self.fitted;self.groups=[np.flatnonzero(z.country_code.eq(g)) for g in sorted(z.country_code.unique())];self.codes=sorted(z.country_code.unique())
  self.within_y=self.y-QD@(QD.T@self.y);self.within_R2=1-float(self.e@self.e)/float(self.within_y@self.within_y);self.full_R2=1-float(self.e@self.e)/float(np.sum((self.y-self.y.mean())**2));self.sse=float(self.e@self.e)
  self.VIF=np.diag(np.linalg.inv(Z.T@Z));self.condition=np.linalg.cond(Z);self.cov=None
  assert np.max(np.abs(self.Q.T@self.e))<1e-7
  if inference:self.infer()
 def infer(self):
  Q=self.Q;L=np.zeros((self.N,self.G,self.p));scores=[];expected=np.zeros((self.p,self.p));self.adjust=[]
  for j,ix in enumerate(self.groups):
   Pgg=np.eye(len(ix))-Q[ix]@Q[ix].T;val,U=eigh(Pgg);pos=val>1e-10;A=(U[:,pos]/np.sqrt(val[pos]))@U[:,pos].T;self.adjust.append(A);expected+=self.B[ix].T@A@Pgg@A@self.B[ix]
   C=A@self.B[ix];scores.append(C.T@self.e[ix]);L[ix,j,:]=C;L[:,j,:]-=Q@(Q[ix].T@C)
  scores=np.array(scores);self.cov=scores.T@scores;self.L=L
  # Exact working-model bias-reduction identity on the real design.
  assert np.allclose(expected,self.B.T@self.B,rtol=2e-6,atol=1e-10)
  self.working_model_covariance_error=float(np.max(np.abs(expected-self.B.T@self.B)))
 def contrast(self,c,label):
  c=np.array(c);b=float(c@self.beta);se=np.sqrt(max(float(c@self.cov@c),0));L=np.einsum('ngp,p->ng',self.L,c);K=L.T@L;df=float(np.trace(K)**2/np.sum(K*K));t=b/se;p=2*stats.t.sf(abs(t),df);crit=stats.t.ppf(.975,df)
  # Design information from the coefficient's OLS influence, not realized residuals.
  a=self.B@c;w=np.array([np.sum(a[ix]**2) for ix in self.groups]);w/=w.sum();order=np.argsort(w)[::-1]
  return dict(term=label,estimate=b,SE_CR2=se,Satterthwaite_df=df,t=t,p_value=p,CI95_low=b-crit*se,CI95_high=b+crit*se,effective_design_clusters=1/float(w@w),largest_information_country=self.codes[order[0]],largest_information_share=w[order[0]],top5_information_share=w[order[:5]].sum())
 def unit(self,name):
  c=np.zeros(self.p);c[self.names.index(name)]=1;return c
 def coefficients(self):return [self.contrast(self.unit(n),n) for n in self.names]
 def joint(self,R,label):
  R=np.asarray(R);q=len(R);theta=R@self.beta;S=R@self.cov@R.T;W=float(theta@np.linalg.solve(S,theta));A=self.B@R.T;Omega=A.T@A;v,U=eigh(Omega);T=(U/np.sqrt(v))@U.T;L=np.einsum('ngp,qp,qr->ngr',self.L,R,T);K=[L[:,:,j].T@L[:,:,j] for j in range(q)];sumK=sum(K);total=float(np.sum(sumK*sumK))
  for i in range(q):
   for j in range(q):
    Kij=L[:,:,i].T@L[:,:,j];total+=float(np.sum(Kij*Kij.T))
  eta=q*(q+1)/total;df2=eta-q+1;F=W*df2/(q*eta);p=stats.f.sf(F,q,df2) if df2>0 else np.nan
  if q==1:assert np.allclose(df2,self.contrast(R[0],label)['Satterthwaite_df'])
  return dict(hypothesis=label,restrictions=q,Wald_chi_square=W,HTZ_F=F,HTZ_df1=q,HTZ_df2=df2,HTZ_p=p)
 def bootstrap(self,R,label,draws=9999,seed=20261016):
  R=np.atleast_2d(R);q=len(R);A=self.B@R.T;theta=R@self.beta;fit0=self.fitted-A@np.linalg.solve(A.T@A,theta);u0=self.y-fit0
  assert np.max(np.abs(R@(self.B.T@fit0)))<1e-7
  U=np.zeros((self.N,self.G))
  for j,ix in enumerate(self.groups):U[ix,j]=u0[ix]
  delta=A.T@U;L=np.einsum('ngp,qp->ngq',self.L,R);tensor=np.einsum('ngq,nh->gqh',L,U);observed=float(theta@np.linalg.solve(R@self.cov@R.T,theta))/q;rng=np.random.default_rng(seed);exceed=0;valid=0
  for start in range(0,draws,250):
   b=min(250,draws-start);weights=rng.choice(np.array([-1.,1.]),size=(self.G,b));th=delta@weights;score=np.einsum('gqh,hb->bgq',tensor,weights,optimize=True);vc=np.einsum('bgq,bgr->bqr',score,score);sol=np.linalg.solve(vc,th.T[:,:,None])[:,:,0];stat=np.sum(th.T*sol,axis=1)/q;exceed+=int(np.sum(stat>=observed));valid+=b
  p=(1+exceed)/(1+valid);return dict(hypothesis=label,restrictions=q,draws=draws,valid_draws=valid,seed=seed,weights='Rademacher',null='Restricted OLS null imposed',studentization='Re-estimated CR2 Wald/q each draw',observed_Wald_F=observed,p_wild_restricted=p,Monte_Carlo_SE=np.sqrt(p*(1-p)/(valid+1)))
