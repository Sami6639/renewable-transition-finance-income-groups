from planb_setup import *
import subprocess,re,os
stage=Path(__file__).parent
if stage.resolve()!=SCRIPTS.resolve():shutil.copy2(stage/'planb_estimate.R',SCRIPTS/'planb_estimate.R')
r=Path(__import__('os').environ.get('RSCRIPT','Rscript'))
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
p=subprocess.Popen([str(r),'--vanilla',str(SCRIPTS/'planb_estimate.R'),str(ROOT)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace',env=env)
lines=[]
for line in p.stdout:
 lines.append(line)
 if line.startswith('MODEL_DONE'):print(line.strip(),flush=True)
code=p.wait();s=''.join(lines);(EV/'R_estimation_full_log.txt').write_text(s,encoding='utf-8');assert code==0,s[-3000:]
for name,body in re.findall(r'BEGIN_EXPORT:([^\r\n]+)\r?\n(.*?)\r?\nEND_EXPORT',s,re.S):(EV/name).write_text(body+'\n',encoding='utf-8')
print('R/clubSandwich estimation finished',flush=True)
