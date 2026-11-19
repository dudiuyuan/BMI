"""Build traceable gene matrices from downloaded official SOFT, preserving raw files."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import tarfile
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
import numpy as np
import pandas as pd
from preprocessing import preprocess


def soft_table(text, kind):
    part = text.split(f'!{kind}_table_begin', 1)[1].split(f'!{kind}_table_end', 1)[0].strip()
    return pd.read_csv(io.StringIO(part), sep='\t', dtype=str, keep_default_na=False)


def prepare(archive, root='data/gse3635', references='data/references'):
    root, refs = Path(root), Path(references)
    raw = root/'raw'
    raw.mkdir(parents=True, exist_ok=True)
    if Path(archive).resolve() != (raw/'GSE3635_RAW.tar').resolve():
        shutil.copy2(archive, raw/'GSE3635_RAW.tar')
    # Never extract arbitrary paths from user tar archives.
    raw_audit = []
    with tarfile.open(raw/'GSE3635_RAW.tar') as tf:
        for member in tf.getmembers():
            if not re.fullmatch(r'GSM\d+\.gpr\.gz', member.name) or not member.isfile():
                raise ValueError(f'Unexpected archive member: {member.name}')
            blob = tf.extractfile(member).read()
            (raw/member.name).write_bytes(blob)
            text = gzip.decompress(blob).decode('latin1')
            lines = text.splitlines()
            start = next(i for i, line in enumerate(lines) if line.startswith('"Block"'))
            tab = pd.read_csv(io.StringIO('\n'.join(lines[start:])), sep='\t')
            raw_audit.append(dict(sample=member.name.split('.')[0], spots=len(tab),
                                  negative_flag_spots=int((tab.Flags < 0).sum()),
                                  unique_raw_ids=tab.ID.nunique()))
    pd.DataFrame(raw_audit).to_csv(root/'raw_gpr_qc.csv', index=False)
    family = gzip.open(raw/'GSE3635_family.soft.gz', 'rt', encoding='utf-8').read()
    platform = soft_table((raw/'GPL1914.soft').read_text(encoding='utf-8'), 'platform')
    if platform.ID.duplicated().any():
        raise ValueError('Ambiguous platform ID mapping')
    platform.to_csv(root/'platform_annotation.csv', index=False)
    # Retain original historical systematic ORFs exactly; no speculative alias remapping.
    mapping = platform[['ID', 'ORF', 'Gene', 'SGDID']].rename(columns={'ID': 'probe_id', 'ORF': 'gene_id'})
    mapping['gene_id'] = mapping.gene_id.str.strip().str.upper()
    mapping['valid_gene_id'] = mapping.gene_id.str.fullmatch(r'Y[A-P][LR](?:\d{3}|X\d{2})[CW](?:-[A-Z])?')
    mapping.to_csv(root/'probe_gene_map.csv', index=False)
    series, flagged, quality, sample_metadata = [], [], [], []
    log_errors = []
    for section in family.split('^SAMPLE = ')[1:]:
        sample = section.splitlines()[0].strip()
        title = re.search(r'^!Sample_title = (.+)$', section, re.M)[1]
        minutes = int(re.search(r'time point (\d+) min', title)[1])
        source_date = re.search(r'\d{4}-\d{2}-\d{2}', title)[0]
        table = soft_table(section, 'sample')
        values = pd.to_numeric(table.VALUE, errors='coerce')
        ratio = pd.to_numeric(table.Ratio, errors='coerce')
        flags = table['Flagged?'].str.upper()
        if not flags.isin(['TRUE', 'FALSE']).all():
            raise ValueError('Unknown sample quality flag')
        mask = flags == 'TRUE'
        good = (ratio > 0) & values.notna()
        log_errors.extend(np.abs(values[good]-np.log10(ratio[good])).tolist())
        s = pd.Series(values.to_numpy(), index=table.ID_REF, name=f't{minutes}')
        if s.index.duplicated().any():
            raise ValueError('Duplicate ID_REF within GEO sample')
        series.append(s)
        flagged.append(pd.Series(mask.to_numpy(), index=table.ID_REF, name=f't{minutes}'))
        quality.append(dict(sample=sample, time_minutes=minutes, probe_rows=len(table),
                            missing_values=int(values.isna().sum()), flagged_values=int(mask.sum())))
        sample_metadata.append(dict(sample_id=sample, time_minutes=minutes, title=title,
                                    source_scan_date=source_date, column=f't{minutes}',
                                    source=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={sample}'))
    samples = pd.DataFrame(sample_metadata).sort_values('time_minutes')
    if samples.time_minutes.tolist() != list(range(0, 121, 10)):
        raise ValueError('Unexpected sample schedule')
    samples.to_csv(root/'samples.csv', index=False)
    pd.DataFrame(quality).sort_values('time_minutes').to_csv(root/'sample_qc.csv', index=False)
    values = pd.concat(series, axis=1)[samples.column.tolist()]
    flag = pd.concat(flagged, axis=1)[samples.column.tolist()]
    values.rename_axis('probe_id').reset_index().to_csv(root/'geo_normalized_probe_values.csv', index=False)
    flag.rename_axis('probe_id').reset_index().to_csv(root/'geo_quality_flags.csv', index=False)
    values = values.mask(flag)
    table = mapping.merge(values.rename_axis('probe_id').reset_index(), on='probe_id', how='right', validate='one_to_one')
    table.to_csv(root/'probe_mapping_audit.csv', index=False)
    good_mapping = table.valid_gene_id.fillna(False).astype(bool)
    genes = table.loc[good_mapping, ['gene_id']+samples.column.tolist()]
    genes = genes.groupby('gene_id', as_index=False, sort=True).median(numeric_only=True)
    genes.to_csv(root/'gene_log10ratio_with_missing.csv', index=False)
    clean, qc, audit = preprocess(genes, 'none', .2)
    clean.to_csv(root/'expression_log10ratio.csv', index=False)
    z, _, _ = preprocess(genes, 'zscore', .2)
    z.to_csv(root/'expression_zscore.csv', index=False)
    audit.to_csv(root/'gene_filter_audit.csv', index=False)
    # Compare independent series matrix to sample tables to detect accidental parsing errors.
    matrix_text = gzip.open(raw/'GSE3635_series_matrix.txt.gz', 'rt').read()
    matrix = soft_table(matrix_text, 'series_matrix').set_index('ID_REF').apply(pd.to_numeric, errors='coerce')
    matrix.columns = samples.column.tolist()
    unmasked = pd.concat(series, axis=1)[samples.column.tolist()]
    aligned = matrix.reindex(unmasked.index)
    matrix_error = float(np.nanmax(np.abs(aligned.to_numpy()-unmasked.to_numpy())))
    if matrix_error > 1e-7:
        raise ValueError('Series matrix and sample SOFT disagree')
    with tarfile.open(refs/'yeastCC_1.50.0.tar.gz') as tf:
        for name in ['orf800.txt', 'CellCycle98.txt', 'createYeastCC.R']:
            (refs/name).write_bytes(tf.extractfile('yeastCC/vignettes/'+name).read())
        for name in ['DESCRIPTION', 'man/orf800.Rd', 'man/spYCCmeta.Rd']:
            (refs/Path(name).name).write_bytes(tf.extractfile('yeastCC/'+name).read())
    known = set((refs/'orf800.txt').read_text().split())
    rows = []
    annotation_lines = (refs/'CellCycle98.txt').read_text(encoding='latin1').splitlines()
    for row in csv.reader(annotation_lines, delimiter='\t'):
        if len(row) <= 8:
            continue
        gene = row[6].strip()
        rule = 'ORF column'
        # Original CellCycle98 has a blank ORF cell for YEL077C, but both
        # YPD and SGD name cells explicitly give YEL077C. Require agreement.
        if not gene and row[7].strip() == row[8].strip() and row[7].strip() in known:
            gene = row[7].strip()
            rule = 'blank ORF; identical YPD and SGD systematic IDs, confirmed in orf800'
        if gene in known:
            rows.append(dict(gene_id=gene, phase=row[3].strip(), symbol=row[8].strip(),
                             phase_order=row[4], cluster_order=row[5], mapping_rule=rule))
    ref = pd.DataFrame(rows)
    if len(known) != 800 or len(ref) != 800 or ref.gene_id.duplicated().any() or set(ref.gene_id) != known:
        raise ValueError('Spellman reference failed 800-gene integrity check')
    ref.to_csv(refs/'spellman800.csv', index=False)
    timestamp = datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
    sources = {
        'GSE3635_RAW.tar': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3635/suppl/GSE3635_RAW.tar',
        'GSE3635_series.soft': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3635&targ=self&form=text&view=full',
        'GSE3635_family.soft.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3635/soft/GSE3635_family.soft.gz',
        'GSE3635_series_matrix.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3635/matrix/GSE3635_series_matrix.txt.gz',
        'GPL1914.soft': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1914&targ=self&form=text&view=full',
        'GSM81064.soft': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM81064&targ=self&form=text&view=full',
        'yeastCC_1.50.0.tar.gz': 'https://bioconductor.org/packages/3.22/data/experiment/src/contrib/yeastCC_1.50.0.tar.gz'}
    manifest = []
    for name, url in sources.items():
        path = refs/name if name.startswith('yeastCC') else raw/name
        manifest.append(dict(file=str(path), url=url, bytes=path.stat().st_size,
                             sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                             retrieved_date='user supplied; original download date unknown' if name == 'GSE3635_RAW.tar' else datetime.fromtimestamp(path.stat().st_mtime, ZoneInfo('Asia/Shanghai')).date().isoformat(),
                             processed_at=timestamp))
    (root/'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    meta = dict(accession='GSE3635', platform='GPL1914', organism='Saccharomyces cerevisiae',
                assay='two-colour spotted DNA/cDNA microarray', time_minutes=samples.time_minutes.tolist(),
                biological_replicates=0, reference='Cy3 asynchronous culture; Cy5 synchronized sample',
                value_description='normalized log ratio using Rosetta Resolver',
                inferred_log_base=10, log_base_evidence='VALUE numerically agrees with log10(Ratio) across all samples',
                max_log10_ratio_residual=max(log_errors), series_matrix_max_error=matrix_error,
                note='Raw GenePix Log Ratio is not the same as GEO normalized VALUE; no raw RNA-seq count model used',
                raw_probe_rows=len(values), unmapped_or_control_rows=int((~good_mapping).sum()),
                duplicate_genes_before_aggregation=int(table.loc[good_mapping, 'gene_id'].duplicated().sum()),
                quality_flags_masked=int(flag.to_numpy().sum()), preprocessing=qc,
                scan_date_caveat='Titles indicate first four samples from 2001-05-03, others from 2001-04-11; time and scan-date effects may be confounded',
                reference_size=len(ref), reference_phase_counts=ref.phase.value_counts().to_dict())
    (root/'metadata.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')
    print(json.dumps(meta, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive', default='../GSE3635_RAW.tar')
    p.add_argument('--output', default='data/gse3635')
    p.add_argument('--references', default='data/references')
    a = p.parse_args()
    prepare(a.archive, a.output, a.references)
