"""Verify all delivered rows against raw UniProt data and test offline reproduction."""
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'


def read(name):
    with (OUT / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    manifest = json.loads((OUT / 'manifest.json').read_text(encoding='utf-8'))
    files = manifest['output_sha256']
    before = {name: sha(OUT / name) for name in files}
    assert before == files, 'Delivered CSV integrity failure'
    ids = {r['gene_id'] for r in read('spellman800_input.csv')}
    original = set((ROOT.parent / 'data/references/orf800.txt').read_text().split())
    assert ids == original and len(ids) == 800, 'Input does not equal original orf800 list'
    entries = {}
    for request in manifest['requests']:
        path = OUT / request['file']
        assert sha(path) == request['sha256']
        for entry in json.loads(path.read_bytes())['results']:
            entries[entry['primaryAccession']] = entry
    expected = {}
    for accession, entry in entries.items():
        assert entry['organism']['taxonId'] == 559292
        gene_names = set()
        for group in entry.get('genes', []):
            for field in ('geneName', 'synonyms', 'orderedLocusNames', 'orfNames'):
                value = group.get(field, [])
                gene_names.update(v['value'] for v in ([value] if isinstance(value, dict) else value))
        for crossref in entry.get('uniProtKBCrossReferences', []):
            if crossref['database'] == 'GO':
                term = next(p['value'] for p in crossref['properties'] if p['key'] == 'GoTerm')
                for gene_id in gene_names & ids:
                    expected[(gene_id, accession, crossref['id'])] = term[2:]
    long = read('GO_long.csv')
    actual = {(r['gene_id'], r['entry_id'], r['GO_code']): r['GO_meaning'] for r in long}
    assert len(actual) == len(long) and actual == expected, 'GO row coverage/content mismatch'
    compact = read('GO.csv')
    assert set(compact[0]) == {'GO_id', 'entry_id', 'GO_meaning', 'GO_code'}
    assert len(compact) == len({(r['entry_id'], r['GO_code']) for r in long})
    assert {(r['entry_id'], r['GO_code'], r['GO_meaning']) for r in compact} == {
        (r['entry_id'], r['GO_code'], r['GO_meaning']) for r in long}
    audit = read('gene_audit.csv')
    assert len(audit) == 800 and {r['gene_id'] for r in audit} == ids
    annotated_ids = {r['gene_id'] for r in long}
    assert {r['gene_id'] for r in audit if r['status'] == 'matched_with_GO'} == annotated_ids
    mapping = read('gene_uniprot_mapping.csv')
    for row in audit:
        pairs = [r for r in mapping if r['gene_id'] == row['gene_id']]
        assert int(row['entry_count']) == len(pairs)
        assert int(row['unique_GO_count']) == len({r['GO_code'] for r in long if r['gene_id'] == row['gene_id']})
        if row['status'] == 'unmatched':
            assert not pairs
        elif row['status'] == 'matched_without_GO':
            assert pairs and row['gene_id'] not in annotated_ids
    subprocess.run([sys.executable, str(ROOT / 'fetch_spellman_go.py'), '--offline'], check=True)
    assert before == {name:sha(OUT / name) for name in files}, 'Offline rebuild changed CSV bytes'
    result = dict(status='passed', original_800_ids_identical=True, raw_response_checksums_valid=True,
                  all_GO_rows_equal_raw_UniProt=True, all_800_genes_accounted_for=True,
                  offline_csv_byte_identical=True, stats=manifest['stats'])
    (OUT / 'verification.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print('Verified raw GO rows, all 800 genes, and byte-identical offline reproduction.')


if __name__ == '__main__':
    verify()
