"""Download UniProt GO annotations for Spellman's 800 genes (Python stdlib only)."""
import argparse
import csv
import hashlib
import json
import platform
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parent
ASPECTS = {'P': 'biological_process', 'F': 'molecular_function', 'C': 'cellular_component'}
FIRST = ['GO_id', 'entry_id', 'GO_meaning', 'GO_code']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def write_csv(path, fields, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def names(entry):
    """Exact names only; no fuzzy symbols, inferred aliases, or cross-strain mapping."""
    result = defaultdict(set)
    for gene in entry.get('genes', []):
        for key in ['geneName', 'synonyms', 'orderedLocusNames', 'orfNames']:
            values = gene.get(key, [])
            if isinstance(values, dict):
                values = [values]
            for value in values:
                result[value['value']].add(key)
    return result


def query_url(query):
    params = dict(query=query, format='json', fields='accession,gene_names,organism_id,go',
                  size=500, sort='accession asc')
    return 'https://rest.uniprot.org/uniprotkb/search?' + urllib.parse.urlencode(params)


class Snapshot:
    def __init__(self, raw, offline):
        self.raw, self.offline, self.requests, self.releases = raw, offline, [], set()
        raw.mkdir(parents=True, exist_ok=True)

    def fetch(self, url):
        key = hashlib.sha256(url.encode()).hexdigest()
        payload, metadata = self.raw / (key + '.json'), self.raw / (key + '.meta.json')
        if payload.exists() and metadata.exists():
            meta = json.loads(metadata.read_text(encoding='utf-8'))
            if meta['url'] != url or meta['sha256'] != sha(payload):
                raise ValueError(f'Cache integrity failed: {payload}')
        else:
            if self.offline:
                raise FileNotFoundError(f'Missing cached response for {url}')
            for attempt in range(5):
                try:
                    request = urllib.request.Request(url, headers={'User-Agent': 'Spellman800-GO/1.0 Python-urllib'})
                    with urllib.request.urlopen(request, timeout=120) as response:
                        body = response.read()
                        headers = dict(response.headers)
                    json.loads(body)  # Do not cache an invalid response.
                    break
                except Exception:
                    if attempt == 4:
                        raise
                    time.sleep(2 ** attempt)
            payload.write_bytes(body)
            meta = dict(url=url, retrieved_at_utc=datetime.now(timezone.utc).isoformat(),
                        headers=headers, sha256=sha(payload))
            dump(metadata, meta)
        h = {k.lower(): v for k, v in meta['headers'].items()}
        release = h.get('x-uniprot-release')
        if not release:
            raise ValueError('Response missing UniProt release header')
        self.releases.add(release)
        if len(self.releases) != 1:
            raise ValueError('Mixed UniProt releases: use a new output directory')
        self.requests.append(dict(meta, file=str(payload.relative_to(self.raw.parent))))
        return json.loads(payload.read_bytes()), h

    def search(self, query):
        url, entries, expected = query_url(query), [], None
        visited = set()
        while url:
            if url in visited:
                raise ValueError('Pagination cycle')
            visited.add(url)
            data, headers = self.fetch(url)
            total = int(headers['x-total-results'])
            expected = total if expected is None else expected
            if total != expected:
                raise ValueError('Result count changed during pagination')
            entries.extend(data['results'])
            print(f'UniProt: {len(entries)}/{expected} entries', flush=True)
            match = re.search(r'<([^>]+)>;\s*rel="next"', headers.get('link', ''))
            url = match[1] if match else None
        if len(entries) != expected or len({e['primaryAccession'] for e in entries}) != expected:
            raise ValueError('Incomplete or duplicate pagination')
        return entries


def build(input_path, out, offline):
    out.mkdir(parents=True, exist_ok=True)
    with input_path.open(encoding='utf-8-sig', newline='') as f:
        reference = list(csv.DictReader(f))
    ids = [r['gene_id'] for r in reference]
    if len(ids) != 800 or len(set(ids)) != 800:
        raise ValueError('Expected exactly 800 unique Spellman ORF identifiers')
    saved_input = out / 'spellman800_input.csv'
    if saved_input.exists() and sha(saved_input) != sha(input_path):
        raise ValueError('Input changed: choose a new output directory')
    saved_input.write_bytes(input_path.read_bytes())
    snapshot = Snapshot(out / 'raw', offline)
    primary_query = '(proteome:UP000002311) AND (organism_id:559292)'
    entries = {e['primaryAccession']: e for e in snapshot.search(primary_query)}
    matched_names = set().union(*(set(names(e)) for e in entries.values()))
    missing = sorted(set(ids) - matched_names)
    fallback_queries = []
    for start in range(0, len(missing), 20):
        group = missing[start:start + 20]
        query = '(organism_id:559292) AND (' + ' OR '.join('gene_exact:' + g for g in group) + ')'
        fallback_queries.append(query)
        for entry in snapshot.search(query):
            entries[entry['primaryAccession']] = entry
    index = defaultdict(list)
    for entry in entries.values():
        if entry['organism']['taxonId'] != 559292:
            raise ValueError('Unexpected organism')
        for gene_id in set(names(entry)) & set(ids):
            index[gene_id].append(entry)
    entry_genes = defaultdict(set)
    for gene_id, hits in index.items():
        for entry in hits:
            entry_genes[entry['primaryAccession']].add(gene_id)
    long_rows, mapping, audit = [], [], []
    for ref in reference:
        gene_id = ref['gene_id']
        hits = sorted(index[gene_id], key=lambda e: e['primaryAccession'])
        gene_terms = set()
        for entry in hits:
            accession = entry['primaryAccession']
            gos = {}
            for crossref in entry.get('uniProtKBCrossReferences', []):
                if crossref['database'] != 'GO':
                    continue
                props = {p['key']: p['value'] for p in crossref.get('properties', [])}
                term = props.get('GoTerm', '')
                if not re.fullmatch(r'GO:\d{7}', crossref['id']) or term[:2] not in ['P:', 'F:', 'C:']:
                    raise ValueError(f'Unexpected GO annotation: {crossref}')
                item = gos.setdefault(crossref['id'], dict(term=term, evidence=set(), references=set(), eco=set()))
                if item['term'] != term:
                    raise ValueError('Conflicting GO term names')
                item['evidence'].add(props.get('GoEvidenceType', ''))
                for ev in crossref.get('evidences', []):
                    item['eco'].add(ev.get('evidenceCode', ''))
                    if ev.get('source') and ev.get('id'):
                        item['references'].add(ev['source'] + ':' + ev['id'])
            for code, item in sorted(gos.items()):
                gene_terms.add(code)
                long_rows.append(dict(entry_id=accession, GO_meaning=item['term'][2:], GO_code=code,
                    gene_id=gene_id, symbol_spellman=ref['symbol'], phase_spellman=ref['phase'],
                    GO_aspect=ASPECTS[item['term'][0]],
                    GO_evidence=';'.join(sorted(item['evidence'] - {''})),
                    ECO_codes=';'.join(sorted(item['eco'] - {''})),
                    references=';'.join(sorted(item['references'])),
                    entry_type=entry['entryType'], mapping_fields=';'.join(sorted(names(entry)[gene_id])),
                    spellman_genes_for_entry=';'.join(sorted(entry_genes[accession])),
                    uniprot_url='https://www.uniprot.org/uniprotkb/' + accession,
                    uniprot_release=next(iter(snapshot.releases))))
            mapping.append(dict(gene_id=gene_id, entry_id=accession, entry_type=entry['entryType'],
                mapping_fields=';'.join(sorted(names(entry)[gene_id])), GO_count=len(gos)))
        status = 'unmatched' if not hits else ('matched_with_GO' if gene_terms else 'matched_without_GO')
        audit.append(dict(gene_id=gene_id, symbol_spellman=ref['symbol'], phase_spellman=ref['phase'],
            status=status, entry_count=len(hits), entry_ids=';'.join(e['primaryAccession'] for e in hits),
            unique_GO_count=len(gene_terms), multiple_entries=len(hits) > 1))
        audit[-1]['shared_entry'] = any(len(entry_genes[e['primaryAccession']]) > 1 for e in hits)
    long_rows.sort(key=lambda r: (r['gene_id'], r['entry_id'], r['GO_code']))
    for number, row in enumerate(long_rows, 1):
        row['GO_id'] = number
    keys = [(r['gene_id'], r['entry_id'], r['GO_code']) for r in long_rows]
    if len(keys) != len(set(keys)):
        raise ValueError('Duplicate gene-entry-GO rows')
    fields = FIRST + ['gene_id', 'symbol_spellman', 'phase_spellman', 'GO_aspect', 'GO_evidence',
        'ECO_codes', 'references', 'entry_type', 'mapping_fields', 'spellman_genes_for_entry', 'uniprot_url', 'uniprot_release']
    write_csv(out / 'GO_long.csv', fields, long_rows)
    # Exact four-column layout of the supplied example. GO_id is a row number.
    compact = { (r['entry_id'], r['GO_code']): r for r in long_rows }
    compact_rows = [dict(compact[key], GO_id=i) for i, key in enumerate(sorted(compact), 1)]
    write_csv(out / 'GO.csv', FIRST, compact_rows)
    write_csv(out / 'gene_uniprot_mapping.csv', ['gene_id','entry_id','entry_type','mapping_fields','GO_count'], mapping)
    audit_fields = ['gene_id','symbol_spellman','phase_spellman','status','entry_count','entry_ids','unique_GO_count','multiple_entries','shared_entry']
    write_csv(out / 'gene_audit.csv', audit_fields, audit)
    write_csv(out / 'genes_needing_review.csv', audit_fields,
              [r for r in audit if r['status'] != 'matched_with_GO' or r['multiple_entries'] or r['shared_entry']])
    counts = dict(Counter(r['status'] for r in audit))
    stats = dict(input_genes=800, statuses=counts, multiple_entry_genes=sum(r['multiple_entries'] for r in audit),
                 shared_entry_genes=sum(r['shared_entry'] for r in audit),
                 shared_entries={k: sorted(v) for k, v in sorted(entry_genes.items()) if len(v) > 1},
                 long_rows=len(long_rows), compact_rows=len(compact_rows),
                 unique_GO_terms=len({r['GO_code'] for r in long_rows}),
                 GO_aspect_rows=dict(Counter(r['GO_aspect'] for r in long_rows)))
    files = ['GO_long.csv','GO.csv','gene_uniprot_mapping.csv','gene_audit.csv','genes_needing_review.csv','spellman800_input.csv']
    dump(out / 'manifest.json', dict(stats=stats, python=platform.python_version(),
        script_sha256=sha(Path(__file__)), input_sha256=sha(input_path),
        primary_query=primary_query, fallback_queries=fallback_queries, taxon_id=559292,
        uniprot_release=next(iter(snapshot.releases)), requests=snapshot.requests,
        output_sha256={name:sha(out / name) for name in files},
        policy='Exact gene names in S. cerevisiae S288C only; keep all matches and evidence types; no GO ancestor expansion; no imputation'))
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=BASE.parent / 'data/references/spellman800.csv')
    parser.add_argument('--output', type=Path, default=BASE / 'results')
    parser.add_argument('--offline', action='store_true', help='Rebuild from verified raw response cache without network')
    args = parser.parse_args()
    build(args.input.resolve(), args.output.resolve(), args.offline)
