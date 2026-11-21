"""Derive exposures, covariates and outcomes for both databases.

Reads data/mimic.duckdb and data/eicu.duckdb created by 01_extract_mimic.py and 02_extract_eicu.py.
Writes data/analytic.csv (one row per ICU stay) and data/pairs.csv (simultaneous SpO2/FiO2 and PaO2/FiO2 pairs).

Time zero is ICU admission. Hourly SpO2 (median per hour) is paired with the most recent FiO2 within 4 h;
FiO2 for low-flow oxygen without a charted value is estimated as 21% + 4% per L/min; 21% when no oxygen record.
"""
import duckdb,pandas as pd,numpy as np,re
def N(x):
    x=x.copy();x['id']=x['id'].astype('int64');x['h']=x['h'].astype('float64');return x
def one_fio2(fi):
    # one FiO2 value per patient and time point; the highest value wins when two records share a timestamp
    return N(fi).groupby(['id','h'],as_index=False).fio2.max()


def hourly_sf(sp,fi):
    sp=N(sp);fi=N(fi)
    sp=sp.copy();sp['hb']=np.floor(sp.h);sp=sp.groupby(['id','hb']).spo2.median().reset_index().rename(columns={'hb':'h'})
    sp=sp.sort_values(['h','id'],kind='mergesort');fi=fi.sort_values(['h','id'],kind='mergesort')
    x=pd.merge_asof(sp,fi,on='h',by='id',direction='backward',tolerance=4.0);x['fio2']=x.fio2.fillna(21);x['sf']=x.spo2/(x.fio2/100);return x
def win(x,col,a,b,fn='min'):
    return x[x.h.between(a,b,inclusive='left')].groupby('id')[col].agg(fn)
def pairs(pa,sx):
    pa=N(pa);sx=N(sx)
    s=sx[sx.spo2<=97][['id','h','spo2','fio2']].sort_values(['h','id'],kind='mergesort');p=pa.sort_values(['h','id'],kind='mergesort')
    q=pd.merge_asof(p,s.rename(columns={'fio2':'f2'}),on='h',by='id',direction='nearest',tolerance=1.0).dropna(subset=['spo2'])
    q['sf']=q.spo2/(q.fio2/100);return q[['id','h','pf','sf']]
def finish(c,sx,pa):
    for k,(a,b) in {'1':(0,24),'3':(48,72),'2':(24,48)}.items():
        c['sf'+k]=c.id.map(win(sx,'sf',a,b));c['pf'+k]=c.id.map(win(pa,'pf',a,b))
    c['sf1_97']=c.id.map(win(sx[sx.spo2<=97],'sf',0,24));c['sf3_97']=c.id.map(win(sx[sx.spo2<=97],'sf',48,72))
    c['spo2min']=c.id.map(win(sx,'spo2',0,24))
    return c
def M():
    c0=duckdb.connect('data/mimic.duckdb',read_only=True);Q=lambda s:c0.execute(s).df()
    c=Q("select * from co");c['id']=c.stay_id;c['T']=pd.to_datetime(c.intime)
    c=c[c.age>=18].copy()
    ce=Q("select * from ce").rename(columns={'stay_id':'id'}).merge(c[['id','T']]);ce['h']=(pd.to_datetime(ce.charttime)-ce['T']).dt.total_seconds()/3600
    sp=ce[(ce.itemid==220277)&ce.valuenum.between(60,100)][['id','h']].assign(spo2=ce.valuenum)
    f1=ce[ce.itemid==223835][['id','h','valuenum']].copy();f1['fio2']=np.where(f1.valuenum<=1,f1.valuenum*100,f1.valuenum)
    f2=ce[(ce.itemid==223834)&ce.valuenum.between(0,15)][['id','h','valuenum']].copy();f2['fio2']=21+4*f2.valuenum
    fi=pd.concat([f1[['id','h','fio2']],f2[['id','h','fio2']]]);fi=one_fio2(fi[fi.fio2.between(21,100)])
    sx=hourly_sf(sp,fi)
    lab=Q("select * from lab").merge(c[['hadm_id','id','T']]);lab['h']=(pd.to_datetime(lab.charttime)-lab['T']).dt.total_seconds()/3600
    art=set(lab.specimen_id[(lab.itemid==52033)&(lab.value=='ART.')])
    pa=lab[(lab.itemid==50821)&lab.specimen_id.isin(art)&lab.valuenum.between(20,700)][['id','h','valuenum']].rename(columns={'valuenum':'pao2'}).sort_values(['h','id'],kind='mergesort')
    pa=pd.merge_asof(N(pa).sort_values(['h','id'],kind='mergesort'),N(fi).sort_values(['h','id'],kind='mergesort'),on='h',by='id',direction='backward',tolerance=4.0);pa['fio2']=pa.fio2.fillna(21);pa['pf']=pa.pao2/(pa.fio2/100)
    pa=N(pa);c['id']=c['id'].astype('int64');c=finish(c,sx,pa);pr=pairs(pa,sx)
    def lv(it,fn,rng,a=-6,b=24):
        d=lab[(lab.itemid==it)&lab.h.between(a,b)&lab.valuenum.between(*rng)];return c.id.map(d.groupby('id').valuenum.agg(fn))
    c['lactate']=lv(50813,'max',(0.1,30));c['creat']=lv(50912,'max',(0.1,20));c['bun']=lv(51006,'max',(1,250));c['wbc']=lv(51301,'max',(0.1,200))
    c['hb']=lv(51222,'min',(3,20));c['plt']=lv(51265,'min',(1,1500));c['bili']=lv(50885,'max',(0.1,50));c['alb']=lv(50862,'min',(1,6));c['ca']=lv(50893,'min',(3,15));c['glu']=lv(50931,'max',(20,2000))
    c['paco2']=c.id.map(lab[(lab.itemid==50818)&lab.specimen_id.isin(art)&lab.h.between(-6,24)].groupby('id').valuenum.max())
    def cv(items,fn,rng):
        d=ce[ce.itemid.isin(items)&ce.h.between(-6,24)&ce.valuenum.between(*rng)];return c.id.map(d.groupby('id').valuenum.agg(fn))
    c['map_min']=cv([220052,220181],'min',(20,200));c['hr']=cv([220045],'max',(20,250));c['rr']=cv([220210,224690],'max',(4,70))
    g=ce[ce.itemid.isin([220739,223900,223901])&ce.h.between(-6,24)].pivot_table(index=['id','charttime'],columns='itemid',values='valuenum').dropna().sum(axis=1).groupby('id').min();c['gcs']=c.id.map(g)
    w=ce[(ce.itemid==226512)&ce.valuenum.between(30,300)].groupby('id').valuenum.first();c['wt']=c.id.map(w)
    va=Q("select * from inp").rename(columns={'stay_id':'id'}).merge(c[['id','T']]);va['h']=(pd.to_datetime(va.starttime)-va['T']).dt.total_seconds()/3600
    c['vaso']=c.id.isin(va.id[va.h.between(-6,24)]).astype(int)
    iv=Q("select * from imv").rename(columns={'stay_id':'id'}).merge(c[['id','T']]);iv['s']=(pd.to_datetime(iv.vs)-iv['T']).dt.total_seconds()/3600;iv['e']=(pd.to_datetime(iv.ve)-iv['T']).dt.total_seconds()/3600
    c['imv24']=c.id.isin(iv.id[(iv.s<=24)&(iv.e>=0)]).astype(int);c['imv0']=c.id.isin(iv.id[(iv.s<=6)&(iv.e>=0)]).astype(int)
    c['imv_after72']=c.id.isin(iv.id[(iv.s>72)]).astype(int);c['imv_0_72']=c.id.isin(iv.id[(iv.s<=72)&(iv.e>=0)]).astype(int);c['imv_new72']=c.id.isin(iv.id[(iv.s>72)&(iv.s<=672)]).astype(int)
    rr=Q("select * from rrt").rename(columns={'stay_id':'id'}).merge(c[['id','T']]);rr['h']=(pd.to_datetime(rr.starttime)-rr['T']).dt.total_seconds()/3600
    c['rrt7']=c.id.isin(rr.id[rr.h.between(0,168)]).astype(int);c['rrt_pre72']=c.id.isin(rr.id[rr.h<=72]).astype(int);c['rrt_new']=c.id.isin(rr.id[(rr.h>72)&(rr.h<=168)]).astype(int)
    dx=Q("select * from dx")
    def has(a,b): m=((dx.icd_version==10)&dx.icd_code.str.match(a))|((dx.icd_version==9)&dx.icd_code.str.match(b));return c.hadm_id.isin(dx.hadm_id[m]).astype(int)
    D=dict(copd=('J4[34]','49[126]'),chf=('I50','428'),ckd=('N18','585'),liver=('K70[34]|K74|K76[67]','571[256]'),dm=('E1[0-4]','250'),cancer=('C[0-8]|C9[0-7]','1[4-9][0-9]|20[0-8]'),
           alc=('K852|F10|K70','3030|3039|305[0]|571[0-3]'),bil=('K851|K80|K83','574|5761'),hld=('E78[1-5]','272[0-4]'),ards=('J80','51882'))
    for k,v in D.items(): c[k]=has(*v)
    c['male']=(c.gender=='M').astype(int);c['sq1']=(c.sq==1).astype(int)
    dt=pd.to_datetime(c.deathtime).fillna(pd.to_datetime(c.dod)+pd.Timedelta(hours=12));c['tdeath']=(dt-c['T']).dt.total_seconds()/86400
    c['dhosp']=c.hospital_expire_flag.astype(int);c['icu_h']=c.los*24;c['hosp_los']=(pd.to_datetime(c.dischtime)-c['T']).dt.total_seconds()/86400
    c['db']='MIMIC-IV';c['has90']=1;c['hosp']=0
    return c,pr
def E():
    c0=duckdb.connect('data/eicu.duckdb',read_only=True);Q=lambda s:c0.execute(s).df()
    c=Q("select * from co");c['id']=c.s;c=c[c.age>=18].copy()
    vp=Q("select s id,o,sao2 from vp where sao2 between 60 and 100");vp['h']=vp.o/60;sp=vp[['id','h']].assign(spo2=vp.sao2)
    rc=Q("select s id,o,k,v from rc");rc['h']=rc.o/60;rc['x']=pd.to_numeric(rc.v.astype(str).str.extract(r'([\d\.]+)')[0],errors='coerce')
    a=rc[rc.k.isin(['FiO2','FIO2 (%)','Set Fraction of Inspired Oxygen (FIO2)','O2 Percentage'])].copy();a['fio2']=np.where(a.x<=1,a.x*100,a.x)
    b=rc[rc.k.isin(['LPM O2','Oxygen Flow Rate'])&rc.x.between(0,15)].copy();b['fio2']=21+4*b.x
    lab=Q("select s id,o,ln,v from lab");lab['h']=lab.o/60;lab['v']=pd.to_numeric(lab.v,errors='coerce')
    l2=lab[lab.ln=='FiO2'].copy();l2['fio2']=np.where(l2.v<=1,l2.v*100,l2.v)
    fi=pd.concat([a[['id','h','fio2']],b[['id','h','fio2']],l2[['id','h','fio2']]]);fi=one_fio2(fi[fi.fio2.between(21,100)])
    sx=hourly_sf(sp,fi)
    pa=lab[(lab.ln=='paO2')&lab.v.between(20,700)][['id','h','v']].rename(columns={'v':'pao2'}).sort_values(['h','id'],kind='mergesort')
    pa=pd.merge_asof(N(pa).sort_values(['h','id'],kind='mergesort'),N(fi).sort_values(['h','id'],kind='mergesort'),on='h',by='id',direction='backward',tolerance=4.0);pa['fio2']=pa.fio2.fillna(21);pa['pf']=pa.pao2/(pa.fio2/100)
    pa=N(pa);c['id']=c['id'].astype('int64');c=finish(c,sx,pa);pr=pairs(pa,sx)
    def lv(n,fn,rng):
        d=lab[(lab.ln==n)&lab.h.between(-6,24)&lab.v.between(*rng)];return c.id.map(d.groupby('id').v.agg(fn))
    c['lactate']=lv('lactate','max',(0.1,30));c['creat']=lv('creatinine','max',(0.1,20));c['bun']=lv('BUN','max',(1,250));c['wbc']=lv('WBC x 1000','max',(0.1,200))
    c['hb']=lv('Hgb','min',(3,20));c['plt']=lv('platelets x 1000','min',(1,1500));c['bili']=lv('total bilirubin','max',(0.1,50));c['alb']=lv('albumin','min',(1,6));c['ca']=lv('calcium','min',(3,15));c['glu']=lv('glucose','max',(20,2000))
    c['paco2']=lv('paCO2','max',(10,150))
    vp2=Q("select s id,o,mbp,hr,rr from vp");vp2['h']=vp2.o/60;v24=vp2[vp2.h.between(-6,24)]
    va_=Q("select s id,o,nimap from vap where nimap between 20 and 200");va_['h']=va_.o/60;mm_=pd.concat([v24[v24.mbp.between(20,200)][['id','mbp']],va_[va_.h.between(-6,24)].rename(columns={'nimap':'mbp'})[['id','mbp']]])
    c['map_min']=c.id.map(mm_.groupby('id').mbp.min());c['hr']=c.id.map(v24[v24.hr.between(20,250)].groupby('id').hr.max());c['rr']=c.id.map(v24[v24.rr.between(4,70)].groupby('id').rr.max())
    ap=Q("select * from aps");ap['gcs']=np.where(ap.meds==1,np.nan,ap[['eyes','motor','verbal']].where(ap[['eyes','motor','verbal']]>0).sum(axis=1,min_count=3));c['gcs']=c.id.map(ap.groupby('s').gcs.first())
    c['wt']=c.w.where(c.w.between(30,300))
    inf=Q("select s id,o,dn from inf");inf['h']=inf.o/60
    c['vaso']=c.id.isin(inf.id[inf.h.between(-6,24)&inf.dn.str.lower().str.contains('norepi|epineph|dopamine|vasopressin|phenylephrine|levophed',regex=True)]).astype(int)
    ep=Q("select * from ep");ep['s_']=ep.vs/60;ep['e']=ep.ve/60
    c['imv24']=c.id.isin(ep.s[(ep.s_<=24)&(ep.e>=0)]).astype(int);c['imv0']=c.id.isin(ep.s[(ep.s_<=6)&(ep.e>=0)]).astype(int);c['imv_after72']=c.id.isin(ep.s[ep.s_>72]).astype(int);c['imv_0_72']=c.id.isin(ep.s[(ep.s_<=72)&(ep.e>=0)]).astype(int);c['imv_new72']=c.id.isin(ep.s[(ep.s_>72)&(ep.s_<=672)]).astype(int)
    tr=Q("select * from tr");c['rrt7']=c.id.isin(tr.s[tr.o.between(0,10080)]).astype(int);c['rrt_pre72']=c.id.isin(tr.s[tr.o<=4320]).astype(int);c['rrt_new']=c.id.isin(tr.s[(tr.o>4320)&(tr.o<=10080)]).astype(int)
    dx=Q("select s,lower(ds) ds from dx");ph=Q("select s,lower(p) p from ph")
    def has(rx): return c.id.isin(set(dx.s[dx.ds.str.contains(rx,regex=True)])|set(ph.s[ph.p.str.contains(rx,regex=True)])).astype(int)
    D=dict(copd='copd|chronic obstructive',chf='congestive heart failure|chf',ckd='renal insufficiency|chronic kidney|ckd',liver='cirrhosis',dm='diabetes',cancer='cancer|malignan|carcinoma|lymphoma|leukemia',
           alc='alcohol|etoh',bil='gallstone|cholelith|biliary|choledoch',hld='hypertriglycerid|hyperlipid',ards='ards|acute respiratory distress')
    for k,v in D.items(): c[k]=has(v)
    c['male']=(c.gender=='Male').astype(int);c['sq1']=c.adx.fillna('').str.lower().str.contains('pancreatitis').astype(int)
    c['tdeath']=np.where(c.hst=='Expired',c.ho/1440,np.nan);c['dhosp']=(c.hst=='Expired').astype(int);c['icu_h']=c.uo/60;c['hosp_los']=c.ho/1440
    c['db']='eICU-CRD';c['has90']=0;c['hosp']=c.hospitalid;c['tdisch']=c.ho/1440
    return c,pr
cm,pm=M();ce_,pe=E()
cols=['db','id','age','male','copd','chf','ckd','liver','dm','cancer','alc','bil','hld','ards','sq1','imv0','imv24','imv_after72','vaso','map_min','hr','rr','gcs','wt','lactate','creat','bun','wbc','hb','plt','bili','alb','ca','glu','paco2',
 'sf1','sf2','sf3','pf1','pf2','pf3','sf1_97','sf3_97','spo2min','tdeath','dhosp','icu_h','hosp_los','rrt7','rrt_pre72','rrt_new','imv_0_72','imv_new72','has90','hosp']
cm['tdisch']=cm.hosp_los
d=pd.concat([cm[cols+['tdisch']],ce_[cols+['tdisch']]],ignore_index=True)
d.to_csv('data/analytic.csv',index=False)
pairs_all=pd.concat([pm.assign(db='MIMIC-IV'),pe.assign(db='eICU-CRD')])
pairs_all.to_csv('data/pairs.csv',index=False)
print(d.groupby('db').agg(n=('id','size'),sf1=('sf1','count'),pf1=('pf1','count')))
