#!/usr/bin/env python3
"""Summarize raw-read QC from FastQC zip/data and runinfo."""
from __future__ import annotations
import argparse, csv, re, zipfile
from pathlib import Path

def parse_fastqc_zip(zpath: Path):
    out = {}
    with zipfile.ZipFile(zpath) as z:
        # find fastqc_data.txt
        names = [n for n in z.namelist() if n.endswith('fastqc_data.txt')]
        if not names:
            return out
        text = z.read(names[0]).decode('utf-8', errors='replace')
    # Total Sequences
    m = re.search(r'Total Sequences\t(\d+)', text)
    if m: out['total_sequences'] = int(m.group(1))
    m = re.search(r'Sequence length\t([^\n]+)', text)
    if m: out['sequence_length'] = m.group(1).strip()
    m = re.search(r'%GC\t([0-9.]+)', text)
    if m: out['pct_gc'] = float(m.group(1))
    # Per base N content mean
    if '>>Per base N content' in text:
        block = text.split('>>Per base N content')[1].split('>>END_MODULE')[0]
        vals=[]
        for line in block.splitlines():
            if line.startswith('#') or not line.strip() or line.startswith('>>'): continue
            parts=line.split('\t')
            if len(parts)>=2:
                try: vals.append(float(parts[1]))
                except: pass
        if vals: out['mean_pct_N'] = sum(vals)/len(vals)
    # Adapter content flag
    m = re.search(r'>>Adapter Content\t(pass|fail|warn)', text)
    if m: out['adapter_content'] = m.group(1)
    m = re.search(r'>>Per base sequence quality\t(pass|fail|warn)', text)
    if m: out['per_base_quality'] = m.group(1)
    m = re.search(r'>>Per sequence quality scores\t(pass|fail|warn)', text)
    if m: out['per_sequence_quality'] = m.group(1)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--fastqc-dir', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    meta = list(csv.DictReader(open(args.metadata), delimiter='\t'))
    fdir = Path(args.fastqc_dir)
    rows=[]
    for r in meta:
        run=r['Run']
        z1=list(fdir.glob(f'{run}_1_fastqc.zip'))
        z2=list(fdir.glob(f'{run}_2_fastqc.zip'))
        if not z1 or not z2:
            rows.append({
                'sample_id': r['sample_id'], 'Run': run,
                'mets_status': r['mets_status'], 'clinical_group': r['clinical_group'],
                'reads_r1': '', 'reads_r2': '', 'length_r1': '', 'length_r2': '',
                'pct_gc_r1': '', 'pct_gc_r2': '', 'mean_pct_N_r1': '', 'mean_pct_N_r2': '',
                'adapter_content_r1': '', 'adapter_content_r2': '',
                'per_base_quality_r1': '', 'per_base_quality_r2': '',
                'qc_flags': 'fastqc_not_run', 'removed': 'no', 'removal_reason': ''
            })
            continue
        d1=parse_fastqc_zip(z1[0])
        d2=parse_fastqc_zip(z2[0])
        # Soft QC thresholds (document only; do not auto-remove)
        flags=[]
        n=d1.get('total_sequences',0)
        if n < 10000: flags.append('low_depth_warn_<10k')
        if d1.get('adapter_content')=='fail' or d2.get('adapter_content')=='fail':
            flags.append('adapter_content_fail')
        if d1.get('per_base_quality')=='fail' or d2.get('per_base_quality')=='fail':
            flags.append('per_base_quality_fail')
        rows.append({
            'sample_id': r['sample_id'],
            'Run': run,
            'mets_status': r['mets_status'],
            'clinical_group': r['clinical_group'],
            'reads_r1': d1.get('total_sequences',''),
            'reads_r2': d2.get('total_sequences',''),
            'length_r1': d1.get('sequence_length',''),
            'length_r2': d2.get('sequence_length',''),
            'pct_gc_r1': d1.get('pct_gc',''),
            'pct_gc_r2': d2.get('pct_gc',''),
            'mean_pct_N_r1': d1.get('mean_pct_N',''),
            'mean_pct_N_r2': d2.get('mean_pct_N',''),
            'adapter_content_r1': d1.get('adapter_content',''),
            'adapter_content_r2': d2.get('adapter_content',''),
            'per_base_quality_r1': d1.get('per_base_quality',''),
            'per_base_quality_r2': d2.get('per_base_quality',''),
            'qc_flags': ';'.join(flags) if flags else 'ok',
            'removed': 'no',
            'removal_reason': ''
        })
    if not rows:
        raise SystemExit('No samples in metadata')
    fields=list(rows[0].keys())
    with open(args.out,'w',newline='') as f:
        w=csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)
    print(f'Wrote {args.out} n={len(rows)} with_fastqc={sum(1 for r in rows if r["qc_flags"]!="fastqc_not_run")}')


if __name__ == '__main__':
    main()
