from qca_setup import *
import subprocess,re
stage=Path(__file__).parent;shutil.copy2(stage/'qca_precision.R',SCRIPTS/'qca_precision.R')
p=subprocess.run([str(Path(__import__('os').environ.get('RSCRIPT','Rscript'))),'--vanilla',str(SCRIPTS/'qca_precision.R'),str(ROOT)],capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120)
(EV/'R_precision_stdout.txt').write_text(p.stdout,encoding='utf-8');(EV/'R_precision_stderr.txt').write_text(p.stderr,encoding='utf-8');assert p.returncode==0
for name,body in re.findall(r'BEGIN_EXPORT:([^\r\n]+)\r?\n(.*?)\r?\nEND_EXPORT',p.stdout,re.S):(EV/name).write_text(body+'\n',encoding='utf-8')
print('Unmasked full-precision QCA C truth-table fit exported; original display tables preserved')
