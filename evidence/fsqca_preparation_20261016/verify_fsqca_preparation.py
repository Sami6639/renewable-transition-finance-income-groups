from qca_replication_setup import *
import pandas as pd,numpy as np,openpyxl,math
from PIL import Image,ImageDraw
P=json.loads((EV/'proposed_protocol_frozen.json').read_text(encoding='utf-8'));assert sha(EV/'proposed_protocol_frozen.json')==(EV/'protocol_freeze_sha256.txt').read_text()
rep=json.loads((EV/'replication_validation.json').read_text());assert rep['all_tolerances_passed'] and rep['coefficients_compared']==25 and rep['joint_tests_compared']==2
payload=json.loads((EV/'workbook_payload.json').read_text(encoding='utf-8'));count=0
for file,sheets in payload.items():
 wb=openpyxl.load_workbook(OUT/file,read_only=True,data_only=True);assert wb.sheetnames==[s['name'] for s in sheets]
 for s in sheets:
  values=list(wb[s['name']].iter_rows(min_row=6,max_row=5+len(s['rows']),max_col=len(s['columns']),values_only=True))
  for i,r in enumerate(s['rows']):
   for j,c in enumerate(s['columns']):
    a=values[i][j];e=r.get(c);count+=1
    if e is None or e=='':assert a is None or a=='',(s['name'],i,c,a,e)
    elif isinstance(e,(int,float)) and not isinstance(e,bool):assert a is not None and math.isclose(float(a),e,rel_tol=1e-9,abs_tol=1e-8),(s['name'],i,c,a,e)
    else:assert a==e,(s['name'],i,c,a,e)
 wb.close()
c=pd.read_csv(OUT/'fsqca_case_dataset.csv');allcases=pd.read_csv(EV/'all_economy_case_eligibility.csv');assert len(c)==111 and c.country_code.nunique()==111 and not c.country_code.duplicated().any();assert len(allcases)==217 and allcases.eligible.sum()==111
assert c.income_group_2018.value_counts().to_dict()=={'Lower middle income':41,'Upper middle income':35,'Low income':24,'High income':11}
src=ROOT/'reports/final_model_selection_20261016/four_income_group_final_dataset.csv';d=pd.read_csv(src).set_index(['country_code','year'])
for r in c.itertuples():
 z=d.loc[r.country_code];assert z.loc[2018,'group']==r.income_group_2018
 assert math.isclose(r.bank_depth_mean_2018_2020,z.loc[2018:2020,'FD.AST.PRVT.GD.ZS'].mean(),rel_tol=1e-12)
 assert math.isclose(r.RQ_mean_2018_2020,z.loc[2018:2020,'regulatory_quality_estimate'].mean(),rel_tol=1e-12,abs_tol=1e-12)
 assert math.isclose(r.real_GDP_pc_mean_2018_2020,z.loc[2018:2020,'gdp_per_capita_constant_2015_usd'].mean(),rel_tol=1e-12)
 y=z.loc[2021:2023,'irena2025_pdf'];pop=z.loc[2021:2023,'population'];assert y.notna().all() and pop.gt(0).all();assert math.isclose(r.financing_mean_musd2022_2021_2023,y.mean(),rel_tol=1e-12,abs_tol=1e-12);assert math.isclose(r.financing_mean_per_capita_2021_2023,(y*1e6/pop).mean(),rel_tol=1e-12)
 for yr in [2021,2022,2023]:assert math.isclose(getattr(r,'financing_musd2022_'+str(yr)),y.loc[yr],rel_tol=1e-12,abs_tol=1e-12)
assert not any('membership' in col.lower() or 'solution' in col.lower() for col in c.columns)
rec=json.loads((EV/'case_preparation_receipt.json').read_text());assert rec['source_sha256']==sha(src) and rec['protocol_sha256']==sha(EV/'proposed_protocol_frozen.json');assert rec['prepared_utc']>P['frozen_utc']
for r in json.loads((EV/'protected_existing_files.json').read_text()):assert sha(Path(r['path']))==r['sha256']
old=json.loads((ROOT/'reports/irena2025_api_evidence/preexisting_files_manifest.json').read_text(encoding='utf-8'))
for r in old:assert sha(Path(r['path']))==r['sha256']
cross=pd.read_csv(EV/'univariate_crossover_diagnostics.csv');assert (cross.raw_below_crossover+cross.raw_equal_crossover+cross.raw_above_crossover==cross.cases).all()
images=list(EV.glob('preview_*.png'));assert len(images)==13;canvas=Image.new('RGB',(1600,math.ceil(len(images)/2)*225),'#E8ECF1');draw=ImageDraw.Draw(canvas)
for i,p in enumerate(images):
 im=Image.open(p).convert('RGB');im.thumbnail((790,198));x=i%2*800;y=i//2*225;draw.text((x+5,y+3),p.stem,fill='#203047');canvas.paste(im,(x+5,y+23))
canvas.save(EV/'all_sheet_previews.png')
files=['independent_econometric_replication.md','fsqca_calibration_protocol.xlsx','fsqca_case_dataset.csv','fsqca_pre_analysis_report.md']
for name in files:assert (OUT/name).is_file() and (OUT/name).stat().st_size>0
result=dict(status='PASSED',independent_R_implementation='R 4.6.1 / clubSandwich 0.7.0',coefficients_replicated=25,joint_tests_replicated=2,independent_replication_no_material_discrepancy=True,case_rows=111,unique_case_countries=111,baseline_high_income_cases=11,all_217_economy_exclusion_ledger=True,workbook_sheets=13,workbook_cells_verified=count,anchors_frozen_before_new_case_diagnostics=True,case_memberships_calculated=False,truth_table_created=False,fsQCA_minimization_run=False,original_result_files_unchanged=True,protected_prior_files_unchanged=len(old),outputs={n:sha(OUT/n) for n in files})
(EV/'final_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
