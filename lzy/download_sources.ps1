# Official public sources. Run from the project11 directory. Local files are retained.
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force data/gse3635/raw,data/references | Out-Null
$sources = @{
 'data/gse3635/raw/GSE3635_series.soft' = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3635&targ=self&form=text&view=full'
 'data/gse3635/raw/GSE3635_family.soft.gz' = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3635/soft/GSE3635_family.soft.gz'
 'data/gse3635/raw/GSE3635_series_matrix.txt.gz' = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3635/matrix/GSE3635_series_matrix.txt.gz'
 'data/gse3635/raw/GPL1914.soft' = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1914&targ=self&form=text&view=full'
 'data/gse3635/raw/GSM81064.soft' = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM81064&targ=self&form=text&view=full'
 'data/references/yeastCC_1.50.0.tar.gz' = 'https://bioconductor.org/packages/3.22/data/experiment/src/contrib/yeastCC_1.50.0.tar.gz'
}
$events = @()
foreach ($entry in $sources.GetEnumerator()) {
 if (!(Test-Path -LiteralPath $entry.Key)) {
  Invoke-WebRequest -Uri $entry.Value -OutFile $entry.Key
  $events += @{file=$entry.Key; url=$entry.Value; downloaded_at=(Get-Date -Format o)}
 }
}
if ($events.Count -gt 0) {
 $events | ConvertTo-Json -Depth 4 | Set-Content -Encoding utf8 download_events.json
}
