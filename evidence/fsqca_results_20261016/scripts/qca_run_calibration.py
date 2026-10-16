from qca_setup import *
import subprocess,re,sys,numpy as np,pandas as pd
from scipy.special import expit
r=Path(__import__('os').environ.get('RSCRIPT','Rscript'))
stage=Path(__file__).parent
if (stage/'qca_calibrate.R').resolve() != (SCRIPTS/'qca_calibrate.R').resolve(): shutil.copy2(stage/'qca_calibrate.R',SCRIPTS/'qca_calibrate.R')
p=subprocess.run([str(r),'--vanilla',str(SCRIPTS/'qca_calibrate.R'),str(ROOT)],capture_output=True,text=True,encoding='utf-8',errors='replace')
(EV/'R_calibration_stdout.txt').write_text(p.stdout,encoding='utf-8');(EV/'R_calibration_stderr.txt').write_text(p.stderr,encoding='utf-8')
print(p.stdout[-2000:]);print(p.stderr);assert p.returncode==0
for name,body in re.findall(r'BEGIN_EXPORT:([^\r\n]+)\r?\n(.*?)\r?\nEND_EXPORT',p.stdout,re.S):
 (OUT/name if name=='fsqca_calibrated_dataset.csv' else EV/name).write_text(body+'\n',encoding='utf-8')
d=pd.read_csv(OUT/'fsqca_calibrated_dataset.csv');src=pd.read_csv(SOURCE)
plan=json.loads((EV/'approved_execution_plan_frozen.json').read_text())
checks={}
for name,(raw,a) in plan['sets'].items():
 x=src[raw].to_numpy();l,c,u=a;v=expit(np.log(19)*(x-c)/np.where(x<c,c-l,u-c));checks[name]=float(np.max(abs(v-d[name].to_numpy())))
assert max(checks.values())<1e-12
assert d.country_code.tolist()==src.country_code.tolist()
assert sorted(d.income_group_2018.value_counts().tolist())==[11,24,35,41]
for prefix,cols,mean in [('FD',['bank_depth_'+str(y) for y in range(2018,2021)],'bank_depth_mean_2018_2020'),('RQ',['RQ_'+str(y) for y in range(2018,2021)],'RQ_mean_2018_2020'),('GDP',['real_GDP_pc_'+str(y) for y in range(2018,2021)],'real_GDP_pc_mean_2018_2020'),('Yabs',['financing_musd2022_'+str(y) for y in range(2021,2024)],'financing_mean_musd2022_2021_2023'),('Ypc',['financing_per_capita_'+str(y) for y in range(2021,2024)],'financing_mean_per_capita_2021_2023')]:
 assert not src[cols].isna().any().any();assert np.allclose(src[cols].mean(axis=1),src[mean],atol=1e-10,rtol=1e-12)
receipt={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':sha(SOURCE),'calibrated_sha256':sha(OUT/'fsqca_calibrated_dataset.csv'),'plan_sha256':sha(EV/'approved_execution_plan_frozen.json'),'cases':len(d),'income_counts':d.income_group_2018.value_counts().to_dict(),'Python_R_max_errors':checks,'freeze_precedes_truth_tables':True}
if (EV/'calibration_frozen_receipt.json').exists():assert json.loads((EV/'calibration_frozen_receipt.json').read_text())['calibrated_sha256']==receipt['calibrated_sha256']
else:(EV/'calibration_frozen_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt,indent=2))
