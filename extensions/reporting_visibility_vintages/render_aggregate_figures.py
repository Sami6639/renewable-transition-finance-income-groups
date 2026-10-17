"""Regenerate the two new manuscript group-coverage figures from public aggregates."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
out=ROOT/'figures'
out.mkdir(exist_ok=True)
def rows(name):
    with (ROOT/'results'/name).open(encoding='utf-8-sig',newline='') as f:
        return {r.get('income',r.get('group')):r for r in csv.DictReader(f)}
groups=['Low income','Lower middle income','Upper middle income','High income']
labels=['Low','Lower-middle','Upper-middle','High']
v=rows('visibility_coverage.csv')
fig,ax=plt.subplots(figsize=(10.5,5))
y=[100*int(v[g]['observed'])/int(v[g]['grid']) for g in groups]
bars=ax.bar(labels,y)
for b,x in zip(bars,y):ax.text(b.get_x()+b.get_width()/2,x+1,f'{x:.1f}%',ha='center')
ax.set_ylim(0,108);ax.set_ylabel('Published financing record present (%)')
ax.set_title('Visibility of IRENA 2025 financing records (2014–2023)')
fig.text(.5,.01,'Reference grid: 2,170 economy-years; absence of a record is not zero national investment.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.055,1,1))
fig.savefig(out/'visibility_chart.png',dpi=220);plt.close(fig)
v=rows('vintage_coverage.csv')
fig,ax=plt.subplots(figsize=(11.5,5))
x=np.arange(len(groups));y25=[100*(int(v[g]['both'])+int(v[g]['2025 only']))/int(v[g]['n']) for g in groups]
y26=[100*(int(v[g]['both'])+int(v[g]['2026 only']))/int(v[g]['n']) for g in groups]
for offset,series,label in [(-.18,y25,'2025 vintage (2022-USD basis)'),(.18,y26,'2026 vintage (2023-USD basis)')]:
    bars=ax.bar(x+offset,series,width=.36,label=label)
    for b,y in zip(bars,series):ax.text(b.get_x()+b.get_width()/2,y+1,f'{y:.1f}%',ha='center',fontsize=8)
ax.set_xticks(x,labels);ax.set_ylim(0,108);ax.set_ylabel('Country-year cells with a reported value (%)')
ax.set_title('Reporting-perimeter sensitivity across IRENA vintages (2015–2023)')
ax.legend()
fig.text(.5,.01,'Counts refer only to record availability; vintages differ in monetary price bases and scope.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.055,1,1))
fig.savefig(out/'vintage_chart.png',dpi=220);plt.close(fig)
print('Generated',out/'visibility_chart.png',out/'vintage_chart.png')
