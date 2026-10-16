from stage2_setup import *
import subprocess
src=Path(__file__).parent/'stage2_independent.R';dst=SCRIPTS/src.name
if src.resolve()!=dst.resolve():shutil.copy2(src,dst)
files=['R_coefficients.csv','R_confidence_intervals.csv','R_joint_tests.csv','R_metadata.csv','R_sample_membership.csv','R_session_info.txt']
if not (EV/'R_stdout.txt').exists():
 r=subprocess.run([str(Path(__import__('os').environ.get('RSCRIPT','Rscript'))),'--vanilla',str(dst)],capture_output=True,timeout=300)
 (EV/'R_stdout.txt').write_bytes(r.stdout);(EV/'R_stderr.txt').write_bytes(r.stderr);(EV/'R_execution.json').write_text(json.dumps(dict(exit_code=r.returncode,stdout_bytes=len(r.stdout),stderr_bytes=len(r.stderr)),indent=2))
 if r.returncode:raise RuntimeError(r.stderr.decode('utf-8',errors='replace'))
text=(EV/'R_stdout.txt').read_bytes().decode('utf-8',errors='replace').replace('\r\n','\n')
for file in files:
 marker='BEGIN_EXPORT:'+file+'\n';start=text.index(marker)+len(marker);end=text.index('\nEND_EXPORT',start);(EV/file).write_text(text[start:end]+'\n',encoding='utf-8');print('Saved',file,flush=True)
