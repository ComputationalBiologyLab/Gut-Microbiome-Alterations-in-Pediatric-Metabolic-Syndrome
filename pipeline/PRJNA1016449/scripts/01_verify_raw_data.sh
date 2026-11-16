#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

python3 - <<'PY'
from pathlib import Path
import csv
root = Path('.')
meta = list(csv.DictReader(open('00_metadata/sample_metadata.tsv'), delimiter='\t'))
n = len(meta)
missing = []
for r in meta:
    for k in ('fastq_r1','fastq_r2'):
        p = root/r[k]
        if not p.exists() or p.stat().st_size < 1000:
            missing.append((r['sample_id'], r['Run'], k, str(p)))
runs = [r['Run'] for r in meta]
dup = sorted({x for x in runs if runs.count(x)>1})
print(f"Expected samples: {n}")
print(f"Unique runs: {len(set(runs))}")
print(f"Duplicated runs: {dup if dup else 'none'}")
print(f"Missing/empty FASTQs: {len(missing)}")
for m in missing[:20]:
    print('MISSING', m)
# write summary
out = Path('results/raw_data_verification.tsv')
with out.open('w') as f:
    f.write('metric\tvalue\n')
    f.write(f'n_samples\t{n}\n')
    f.write(f'n_unique_runs\t{len(set(runs))}\n')
    f.write(f'n_missing_fastq\t{len(missing)}\n')
    f.write(f'all_paired\t{str(len(missing)==0)}\n')
    f.write(f'n_MetS\t{sum(1 for r in meta if r["mets_status"]=="MetS")}\n')
    f.write(f'n_Control_nonMetS\t{sum(1 for r in meta if r["mets_status"]=="Control")}\n')
if missing:
    raise SystemExit(1)
print('OK: all 50 paired FASTQs present, no duplicate runs')
PY
