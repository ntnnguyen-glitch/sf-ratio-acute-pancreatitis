"""Extract the eICU-CRD v2.0 acute pancreatitis cohort and raw measurements.

Set EICU_DIR to the folder that contains the eICU-CRD v2.0 .csv.gz files.
Output: data/eicu.duckdb
"""
import os
import time
import duckdb

EICU_DIR = os.environ["EICU_DIR"].rstrip("/") + "/"
os.makedirs("data", exist_ok=True)
t0 = time.time()


def log(*msg):
    print(round(time.time() - t0), *msg, flush=True)


con = duckdb.connect("data/eicu.duckdb")
X = con.execute
X("PRAGMA threads=8")


def csv(table, varchar=False):
    if varchar:
        return f"read_csv('{EICU_DIR}{table}.csv.gz', all_varchar=true, strict_mode=false, ignore_errors=true)"
    return f"read_csv('{EICU_DIR}{table}.csv.gz', header=true, auto_detect=true, sample_size=-1, ignore_errors=true)"


X(f"create or replace table pt as select * from {csv('patient')}")
X(f"create or replace table dx_all as select patientunitstayid s, lower(diagnosisstring) ds, icd9code ic from {csv('diagnosis')}")
log("base tables")

# Cohort: APACHE admission diagnosis or problem-list diagnosis of acute pancreatitis; first ICU stay per patient
X("""
create or replace table co as
with a as (
  select distinct s from dx_all
  where ds like 'gastrointestinal|pancreatic disease|pancreatitis%' and ds not like '%chronic%'
  union
  select patientunitstayid from pt where lower(apacheadmissiondx) like '%pancreatitis%'),
p as (
  select t.*, row_number() over (partition by t.uniquepid
                                 order by t.hospitaldischargeyear, t.patienthealthsystemstayid, t.unitvisitnumber) rn
  from pt t join a on a.s = t.patientunitstayid)
select patientunitstayid s, patienthealthsystemstayid hs, hospitalid, gender,
       case when cast(age as varchar) like '>%' then 90 else try_cast(age as int) end age,
       ethnicity, apacheadmissiondx adx, admissionweight w, uniquepid,
       hospitaladmitoffset hao, hospitaldischargeoffset ho, hospitaldischargestatus hst,
       unitdischargeoffset uo, unitdischargestatus ust
from p where rn = 1
""")
log("cohort", X("select count(*) from co").fetchone()[0])

X(f"""
create or replace table vp as
select try_cast(patientunitstayid as int) s, try_cast(observationoffset as int) o, try_cast(sao2 as double) sao2,
       try_cast(heartrate as double) hr, try_cast(respiration as double) rr, try_cast(systemicmean as double) mbp
from {csv('vitalPeriodic', True)}
where try_cast(patientunitstayid as int) in (select s from co) and try_cast(observationoffset as int) between -360 and 10080
""")
log("vitalPeriodic")
X(f"""
create or replace table vap as
select try_cast(patientunitstayid as int) s, try_cast(observationoffset as int) o, try_cast(noninvasivemean as double) nimap
from {csv('vitalAperiodic', True)}
where try_cast(patientunitstayid as int) in (select s from co) and try_cast(observationoffset as int) between -360 and 1440
""")
X(f"""
create or replace table rc as
select try_cast(patientunitstayid as int) s, try_cast(respchartoffset as int) o, respchartvaluelabel k, respchartvalue v
from {csv('respiratoryCharting', True)}
where try_cast(patientunitstayid as int) in (select s from co) and try_cast(respchartoffset as int) between -360 and 10080
  and respchartvaluelabel in ('FiO2', 'FIO2 (%)', 'Set Fraction of Inspired Oxygen (FIO2)', 'LPM O2', 'O2 Percentage',
                              'Oxygen Flow Rate', 'O2 Device', 'Oxygen Delivery Method', 'PEEP')
""")
log("respiratoryCharting")
X(f"""
create or replace table lab as
select patientunitstayid s, labresultoffset o, labname ln, labresult v from {csv('lab')}
where patientunitstayid in (select s from co) and labresultoffset between -1440 and 10080
  and labname in ('paO2', 'FiO2', 'paCO2', 'pH', 'lactate', 'creatinine', 'BUN', 'WBC x 1000', 'platelets x 1000',
                  'Hgb', 'total bilirubin', 'bicarbonate')
""")
log("lab")
# Invasive ventilation episodes recorded in respiratoryCare (prior ventilation start and end offsets)
X(f"""
create or replace table ep as
select distinct patientunitstayid s, priorventstartoffset vs, priorventendoffset ve from {csv('respiratoryCare')}
where patientunitstayid in (select s from co) and priorventstartoffset < priorventendoffset
""")
X(f"""
create or replace table inf as
select patientunitstayid s, infusionoffset o, lower(drugname) dn from {csv('infusionDrug')}
where patientunitstayid in (select s from co) and infusionoffset between -360 and 10080
  and regexp_matches(lower(drugname), 'norepi|levophed|epinephrine|dopamine|phenylephrine|vasopressin')
""")
X("create or replace table dx as select * from dx_all where s in (select s from co)")
X(f"create or replace table ph as select patientunitstayid s, lower(pasthistorypath) p from {csv('pastHistory')} where patientunitstayid in (select s from co)")
X(f"create or replace table aps as select patientunitstayid s, eyes, motor, verbal, meds, dialysis from {csv('apacheApsVar')} where patientunitstayid in (select s from co)")
X(f"""
create or replace table tr as
select patientunitstayid s, treatmentoffset o, treatmentstring t from {csv('treatment')}
where patientunitstayid in (select s from co) and regexp_matches(lower(treatmentstring), 'dialysis|hemofiltration|crrt')
""")
X("drop table dx_all")
log("done")
