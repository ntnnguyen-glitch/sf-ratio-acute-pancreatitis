"""Extract the MIMIC-IV v3.1 acute pancreatitis cohort and raw measurements.

Set MIMIC_DIR to the folder that contains the hosp/ and icu/ subfolders of MIMIC-IV v3.1.
Output: data/mimic.duckdb
"""
import os
import time
import duckdb

MIMIC_DIR = os.environ["MIMIC_DIR"].rstrip("/") + "/"
os.makedirs("data", exist_ok=True)
t0 = time.time()


def log(*msg):
    print(round(time.time() - t0), *msg, flush=True)


con = duckdb.connect("data/mimic.duckdb")
X = con.execute
X("PRAGMA threads=8")


def csv(path):
    return f"read_csv('{MIMIC_DIR}{path}.csv.gz', header=true, auto_detect=true, sample_size=200000, ignore_errors=true)"


# Base tables
X(f"create or replace table icu as select subject_id, hadm_id, stay_id, first_careunit, intime, outtime, los from {csv('icu/icustays')}")
X(f"create or replace table pat as select subject_id, gender, anchor_age, anchor_year, dod from {csv('hosp/patients')}")
X(f"create or replace table adm as select subject_id, hadm_id, admittime, dischtime, deathtime, race, hospital_expire_flag from {csv('hosp/admissions')}")
X(f"create or replace table dx_all as select subject_id, hadm_id, seq_num, cast(icd_code as varchar) icd_code, icd_version from {csv('hosp/diagnoses_icd')}")
X(f"create or replace table pe_all as select stay_id, itemid, starttime, endtime from {csv('icu/procedureevents')} "
  "where itemid in (225792, 225802, 225803, 225805, 224270, 225809, 225955, 225441)")
log("base tables")

# Cohort: acute pancreatitis among the first three diagnoses, first ICU stay per patient
X("""
create or replace table co as
with d as (
  select hadm_id, min(seq_num) sq from dx_all
  where (icd_version = 10 and icd_code like 'K85%') or (icd_version = 9 and icd_code = '5770')
  group by 1),
i as (
  select i.*, d.sq, row_number() over (partition by i.subject_id order by i.intime) rn
  from icu i join d using (hadm_id))
select i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime, i.first_careunit, i.los, i.sq,
       p.gender, p.anchor_age + (year(i.intime) - p.anchor_year) age, p.dod,
       a.deathtime, a.dischtime, a.admittime, a.race, a.hospital_expire_flag
from i join pat p using (subject_id) join adm a using (hadm_id)
where rn = 1 and sq <= 3
""")
log("cohort", X("select count(*) from co").fetchone()[0])

X(f"""
create or replace table ce as
select c.stay_id, c.itemid, c.charttime, c.valuenum, c.value from {csv('icu/chartevents')} c join co m using (stay_id)
where c.itemid in (220277, 223835, 223834, 226732, 220045, 220210, 224690, 220052, 220181, 220739, 223900, 223901,
                   223762, 223761, 226512, 220339, 224696)
  and c.charttime between m.intime - interval 6 hour and m.intime + interval 7 day
""")
log("chartevents")
X(f"""
create or replace table lab as
select l.hadm_id, l.specimen_id, l.itemid, l.charttime, l.valuenum, l.value from {csv('hosp/labevents')} l join co m using (hadm_id)
where l.itemid in (50821, 50818, 50820, 52033, 50813, 50912, 51006, 51301, 51222, 51265, 50885, 50862, 50893, 50931,
                   50956, 50882, 50983, 50971, 50970)
  and l.charttime between m.intime - interval 24 hour and m.intime + interval 7 day
""")
log("labevents")
X(f"""
create or replace table inp as
select i.stay_id, i.itemid, i.starttime, i.endtime, i.amount from {csv('icu/inputevents')} i join co m using (stay_id)
where i.itemid in (221906, 221289, 221662, 222315, 221749, 229617)
  and i.starttime between m.intime - interval 6 hour and m.intime + interval 7 day
""")
# Invasive ventilation episodes (procedureevents 225792) and renal replacement therapy
X("create or replace table imv as select p.stay_id, p.starttime vs, p.endtime ve from pe_all p join co using (stay_id) where p.itemid = 225792")
X("create or replace table rrt as select p.stay_id, p.starttime from pe_all p join co using (stay_id) where p.itemid <> 225792")
X("create or replace table dx as select d.hadm_id, d.icd_code, d.icd_version from dx_all d join co using (hadm_id)")
for tbl in ["dx_all", "pe_all"]:
    X(f"drop table {tbl}")
log("done")
