from pathlib import Path
import subprocess,json
OUT=Path('[PROJECT_ROOT]/fsqca_preparation_20261016');EV=OUT/'audit_evidence'
files=['R_independent_sample.csv','R_clubSandwich_coefficients.csv','R_clubSandwich_confidence_intervals.csv','R_clubSandwich_joint_tests.csv','R_model_metadata.csv','R_session_info.txt','clubSandwich_version.txt']
code="source('"+str(OUT/'replicate_clubsandwich.R').replace('\\','/')+"',encoding='UTF-8'); for(n in c("+','.join("'"+n+"'" for n in files)+")){cat('\\nBEGIN_EXPORT:',n,'\\n',sep='');cat(paste(readLines(file.path(ev,n),warn=FALSE),collapse='\\n'));cat('\\nEND_EXPORT\\n')}"
if not (EV/'R_captured_stdout.txt').exists():
 r=subprocess.run([str(OUT/'tools/R/bin/Rscript.exe'),'--vanilla','-e',code],capture_output=True,timeout=180)
 (EV/'R_captured_stdout.txt').write_bytes(r.stdout);(EV/'R_captured_stderr.txt').write_bytes(r.stderr);print('R exit',r.returncode,'stdout bytes',len(r.stdout),'stderr bytes',len(r.stderr),flush=True)
 if r.returncode:raise RuntimeError(r.stderr.decode('utf-8',errors='replace'))
text=(EV/'R_captured_stdout.txt').read_bytes().decode('utf-8',errors='replace').replace('\r\n','\n')
for n in files:
 marker='BEGIN_EXPORT:'+n+'\n';start=text.index(marker)+len(marker);end=text.index('\nEND_EXPORT',start);(EV/n).write_text(text[start:end]+'\n',encoding='utf-8');print('Saved',n,flush=True)
