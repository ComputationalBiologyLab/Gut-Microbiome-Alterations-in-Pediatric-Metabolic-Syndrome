#!/usr/bin/env python3
"""Summarize raw-read QC from FastQC zip/data (single-end)."""
from __future__ import annotations
import argparse, csv, re, zipfile
from pathlib import Path

def parse_fastqc_zip(zpath: Path):
    out = {}
    with zipfile.ZipFile(zpath) as z:
        names = [n for n in z.namelist() if n.endswith('fastqc_data.txt')]
        if not names:
            return out
        text = z.read(names[0]).decode('utf-8', errors='replace')
    m = re.search(r'Total Sequences\t(\d+)', text)
    if m: out['total_sequences'] = int(m.group(1))
    m = re.search(r'Sequence length\t([^\n]+)', text)
    if m: out['sequence_length'] = m.group(1).strip()
    m = re.search(r'%GC\t([0-9.]+)', text)
    if m: out['pct_gc'] = float(m.group(1))
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
        run = r.get('run_accession') or r.get('Run')
        z1 = list(fdir.glob(f'{run}_1_fastqc.zip')) or list(fdir.glob(f'{run}_fastqc.zip'))
        if not z1:
            rows.append({
                'sample_id': r['sample_id'], 'Run': run,
                'mets_status': r['mets_status'],
                'reads': '', 'length': '', 'pct_gc': '', 'mean_pct_N': '',
                'adapter_content': '', 'per_base_quality': '',
                'qc_flags': 'fastqc_not_run', 'removed': 'no', 'removal_reason': ''
            })
            continue
        d1 = parse_fastqc_zip(z1[0])
        flags=[]
        n = d1.get('total_sequences', 0)
        if n and n < 10000: flags.append('low_depth_warn_<10k')
        if d1.get('adapter_content') == 'fail': flags.append('adapter_content_fail')
        if d1.get('per_base_quality') == 'fail': flags.append('per_base_quality_fail')
        rows.append({
            'sample_id': r['sample_id'],
            'Run': run,
            'mets_status': r['mets_status'],
            'reads': d1.get('total_sequences', ''),
            'length': d1.get('sequence_length', ''),
            'pct_gc': d1.get('pct_gc', ''),
            'mean_pct_N': d1.get('mean_pct_N', ''),
            'adapter_content': d1.get('adapter_content', ''),
            'per_base_quality': d1.get('per_base_quality', ''),
            'qc_flags': ';'.join(flags) if flags else 'ok',
            'removed': 'no',
            'removal_reason': ''
        })
    fields = list(rows[0].keys())
    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)
    print(f'Wrote {args.out} n={len(rows)} with_fastqc={sum(1 for r in rows if r["qc_flags"]!="fastqc_not_run")}')

if __name__ == '__main__':
    main()
