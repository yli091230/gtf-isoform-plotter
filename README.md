# GTF Isoform Plotter

> **Notice:** This software was generated with assistance from OpenAI Codex.
> Review the source code and validate every plot against the original
> annotation data before publication or downstream analysis. The software is
> provided without warranty.

Create publication-ready transcript diagrams from plain or gzip-compressed GTF
files. Output is an editable PDF: labels remain text, while exons, introns,
strand chevrons, axes, highlights, and scale bars remain vector objects.

## 1. Quick start

Python 3.9 or newer is required.

```bash
cd gtf-isoform-plotter
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

Plot every annotated isoform of one gene:

```bash
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type isoforms \
  --genes SEPTIN3 \
  --unit kb \
  --output SEPTIN3.pdf
```

The default plotting mode is `isoforms`. The explicit `--plot-type isoforms`
shown throughout this guide makes the intended workflow clear.

## 2. Features

- Reads standard `.gtf` and `.gtf.gz` files.
- Selects genes by `gene_name` or `gene_id`.
- Draws all requested isoforms from top to bottom and labels them with their
  stable `transcript_id`, such as `ENST...`.
- Shows narrow UTR/noncoding exon segments and taller CDS segments.
- Preserves literal genomic widths for isoform-mode exons, CDS segments,
  affected regions, and splice-boundary shades.
- Adds 5% of the plotted gene span to each side of an isoform x-axis.
- Writes a separate, gene-named, editable PDF for every gene in a multi-gene
  request.
- Provides a separate reference-region mode with one representative transcript
  per gene, arranged into forward- and reverse-strand rows.
- Filters reference tracks by gene biotype; protein-coding genes are the
  default.
- Uses uniformly spaced, coordinate-aligned strand chevrons.
- Displays coordinates in kb or Mb.
- Optionally highlights affected exon regions in red and splice boundaries in
  light blue.
- Marks transcript-only records with orange hatched blocks so they cannot be
  mistaken for known single-exon structures.

## 3. Installation

### 3.1 pip

Use the commands in [Quick start](#1-quick-start). For an editable development
installation with test dependencies:

```bash
python -m pip install -e ".[dev]"
```

### 3.2 Conda

```bash
conda env create -f environment.yml
conda activate gtf-isoform-plotter
```

### 3.3 Run without installation

From the repository root, the compatibility launcher accepts the same options:

```bash
python plot_transcript.py \
  --gtf-file annotation.gtf.gz \
  --plot-type isoforms \
  --genes SEPTIN3
```

## 4. Plotting modes

The two workflows are intentionally separate:

| Mode | Purpose | Required selection |
|---|---|---|
| `isoforms` | Plot transcript isoforms for one or more genes | `--genes` or `--gene-file` |
| `reference` | Plot one representative transcript per gene in a region | `--region` |

### 4.1 Isoform mode

#### 4.1.1 One gene

```bash
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type isoforms \
  --genes SEPTIN3 \
  --output SEPTIN3.pdf
```

Gene IDs are also accepted:

```bash
gtf-isoform-plotter \
  --gtf-file annotation.gtf.gz \
  --plot-type isoforms \
  --genes ENSG00000100167.19 \
  --output SEPTIN3_by_gene_id.pdf
```

Editable example: [isoform_example.pdf](docs/examples/isoform_example.pdf)

#### 4.1.2 Multiple genes

Pass comma-separated names or IDs. Each gene is written to a separate PDF:

```bash
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type isoforms \
  --genes SEPTIN3,DGKZ,PAK6 \
  --unit Mb \
  --output-dir selected_genes
```

```text
selected_genes/
├── SEPTIN3.pdf
├── DGKZ.pdf
└── PAK6.pdf
```

For one gene, `--output custom.pdf` sets the exact filename. For multiple genes,
use `--output-dir`; if it is omitted, gene-named files are written beside the
`--output` path or in the current directory.

#### 4.1.3 Gene table

Use `--gene-file` instead of `--genes` to provide genes and optional
gene-specific splice annotations:

```bash
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type isoforms \
  --gene-file genes_and_splicing.tsv \
  --output-dir selected_isoforms
```

The file has one to three **tab-separated** columns:

```text
gene       affected_exon       splicing_boundary
SEPTIN3    41985984,41986112    41986112,41987206
DGKZ       46369940,46370009    46369363,46371384
PAK6
```

| Column | Required | Meaning |
|---|---|---|
| 1 | Yes | Gene name or gene ID |
| 2 | No | Affected-exon coordinate pair |
| 3 | No | Splice-boundary coordinate pair; column 2 must also be present |

A header is optional. Blank lines and lines beginning with `#` are ignored.
Repeat a gene on multiple rows to apply multiple events; the gene is plotted
once, in the order of its first appearance.

A one-column row produces an ordinary plot. A two-column row adds red
affected-region highlighting. A three-column row adds both the red annotation
and the full-height light-blue boundary shade.

Use exactly one of `--genes` and `--gene-file`. If `--gene-file` annotations
and a global `--splice-file` are both supplied, their annotations are combined.

#### 4.1.4 Transcript filtering

The default is `--transcript-filter all`. Other supported filters are:

| Filter | Selection rule |
|---|---|
| `all` | Every transcript |
| `basic` | GTF tag `basic` |
| `mane` | GTF tag `MANE_Select` |
| `canonical` | GTF tag `Ensembl_canonical` |
| `appris_principal` | Tag beginning with `appris_principal` |
| `protein_coding` | `transcript_type "protein_coding"` |
| `coding` | At least one GTF `CDS` record |
| `noncoding` | No GTF `CDS` records |
| `type:VALUE` | Exact custom `transcript_type` |

Combine filters with commas; combinations use **OR** behavior:

```bash
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type isoforms \
  --genes SEPTIN3 \
  --transcript-filter mane,basic \
  --output SEPTIN3_selected.pdf
```

Aliases `mane_select`, `ensembl_canonical`, and `appris` are accepted. The
command stops with an error if no transcript passes the requested filters.
Filter accuracy depends on the tags and CDS records supplied by the GTF.

#### 4.1.5 UTR, CDS, and coordinate scale

Isoform mode draws every complete exon as a narrow block, then overlays its GTF
`CDS` intervals as taller blocks:

```text
5′ UTR       CDS                    3′ UTR
━━━━━━████████████████████████████━━━
```

- Exon plus CDS records produce narrow UTRs and tall CDS segments.
- Noncoding transcripts with exon records remain narrow.
- If CDS records are absent, coding boundaries are not inferred.
- Exon, CDS, highlight, and splice-boundary widths remain proportional to their
  literal genomic lengths; no minimum visible width is imposed.
- The x-axis includes 5% of the gene span as a flank on each side.

#### 4.1.6 Splice highlighting

`--splice-file` is optional. Omit it to produce a plot without red or
light-blue annotations. The file has two tab-separated coordinate-pair columns:

```text
# affected_exon        splicing_boundary
41985984,41986112      41986112,41987206
```

```bash
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type isoforms \
  --genes SEPTIN3 \
  --splice-file SEPTIN3_splice.tsv \
  --splice-tolerance 2 \
  --output SEPTIN3_splice.pdf
```

Affected-region rules:

1. If both input endpoints match a GTF exon within `--splice-tolerance`
   (default ±2 bp), the complete exon is red.
2. Otherwise, an annotation matches when at least 80% of the supplied interval
   overlaps a GTF exon; only the intersection is red.
3. If no known exon matches, nothing is colored.
4. UTR/CDS height remains visible within the red region.
5. Transcript-only fallback spans are not treated as exons.

The second coordinate pair is drawn once as one continuous, full-height
light-blue box across all transcript rows. Its left endpoint is compared with
known exon starts and its right endpoint with known exon ends; an endpoint
snaps to the nearest corresponding boundary only when it is within
`--splice-tolerance`. Otherwise, the input coordinate is retained. Use
`--splice-tolerance 0` to require exact coordinate matches.

### 4.2 Reference-region mode

Reference mode plots genes overlapping a chromosome interval and selects one
representative transcript for each gene. By default, only protein-coding genes
are included:

```bash
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type reference \
  --region chr22:41900000-42100000 \
  --reference-gene-type protein_coding \
  --unit Mb \
  --output chr22_reference_genes.pdf
```

Commas are accepted:

```text
--region chr22:41,900,000-42,100,000
```

Editable example:
[reference_example.pdf](docs/examples/reference_example.pdf)

#### 4.2.1 Gene-type filtering

Use `--reference-gene-type` to select gene biotypes:

| Filter | Selection rule |
|---|---|
| `protein_coding` | Protein-coding genes; this is the default |
| `all` | Every gene type |
| `noncoding` | Every annotated type other than `protein_coding` |
| `pseudogene` | Gene type containing `pseudogene` |
| `miRNA` | miRNA genes |
| `lncRNA` | `lncRNA` or `lincRNA` genes |
| `type:VALUE` | Exact `gene_type` or `gene_biotype` value |

Filters can be combined with commas and use OR behavior:

```bash
# Plot every annotated gene type
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type reference \
  --region chr22:40110000-42420000 \
  --reference-gene-type all \
  --unit Mb \
  --output chr22_all_reference_genes.pdf

# Plot miRNA and lncRNA genes
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type reference \
  --region chr22:40110000-42420000 \
  --reference-gene-type miRNA,lncRNA \
  --unit Mb \
  --output chr22_ncRNA_reference_genes.pdf

# Plot one exact provider-specific biotype
gtf-isoform-plotter \
  --gtf-file gencode.v38.annotation.gtf.gz \
  --plot-type reference \
  --region chr22:40110000-42420000 \
  --reference-gene-type type:processed_pseudogene \
  --unit Mb \
  --output chr22_processed_pseudogenes.pdf
```

The plotter reads `gene_type`, then `gene_biotype`, and finally
`transcript_type` when the preceding attributes are absent. Genes without
usable biotype metadata are included by `all` but not by a specific type
filter.

#### 4.2.2 Layout

- Forward-strand genes occupy the top row.
- Reverse-strand genes occupy the bottom row.
- Gene labels appear above forward models and below reverse models.
- Nearby labels are staggered to reduce collisions.
- A scale bar is added at the upper right in the requested `--unit`.

#### 4.2.3 Representative transcript selection

Candidates are ranked in this order:

1. Curated RefSeq accessions beginning with `NM_` or `NR_`
2. GENCODE `MANE_Select` transcripts
3. `Ensembl_canonical` transcripts
4. Transcripts tagged `basic`
5. Longest available transcript

The plotter selects one complete transcript and draws that transcript's
overlapping exons. It never assembles a model from exons belonging to different
transcripts.

```text
Find genes overlapping --region
        ↓
Collect candidate transcripts
        ↓
Rank RefSeq/MANE/canonical/basic/length
        ↓
Choose one transcript per gene
        ↓
Plot its exons within the visible region
```

Exons crossing a region boundary are clipped at the plot edge. Exons outside
the interval are not shown. When available, full transcript bounds—not only
visible exons—are used for longest-transcript ranking. Selection can change
with the GTF release or annotation provider.

## 5. Command-line reference

| Option | Meaning |
|---|---|
| `--gtf-file PATH` | Required `.gtf` or `.gtf.gz` input |
| `--plot-type {isoforms,reference}` | Plotting workflow; default `isoforms` |
| `--genes NAMES` | Comma-separated gene names/IDs for isoform mode |
| `--gene-file PATH` | Isoform TSV with gene and optional annotation columns |
| `--region REGION` | `chrom:start-end` selection for reference mode |
| `--transcript-filter LIST` | Comma-separated isoform filters; default `all` |
| `--reference-gene-type LIST` | Comma-separated reference gene types; default `protein_coding` |
| `-o, --output PATH` | Single-gene isoform or reference PDF |
| `--output-dir PATH` | Directory for separate gene-named isoform PDFs |
| `--unit {kb,Mb}` | Coordinate unit; default `kb` |
| `--splice-file PATH` | Optional two-column annotation file |
| `--splice-tolerance BP` | Coordinate tolerance; default `2` |
| `-h, --help` | Display command help |

Mode requirements:

- Isoform mode requires exactly one of `--genes` and `--gene-file`.
- Reference mode requires `--region`.
- `--transcript-filter` and `--output-dir` are isoform-mode options.
- `--reference-gene-type` controls reference-mode gene-biotype selection.

## 6. GTF requirements and limitations

The parser uses standard nine-column GTF records and reads `gene_name` or
`gene_id`, `transcript_id`, chromosome, strand, exon coordinates, CDS
coordinates, transcript type, and supported tags.

Full exon–intron structures require `exon` features. A `transcript` record alone
provides only the outer span:

```text
chr1  source  transcript  131125  135623  .  +  .  transcript_id "TX1";
```

It cannot distinguish between structures such as:

```text
131125  █████────████──────█████████████  135623
131125  ███────────────████──────────████  135623
```

A complete model requires exon records:

```text
chr22  source  transcript  100  900  .  +  .  transcript_id "TX1";
chr22  source  exon        100  200  .  +  .  transcript_id "TX1";
chr22  source  exon        400  500  .  +  .  transcript_id "TX1";
chr22  source  exon        700  900  .  +  .  transcript_id "TX1";
```

To expose missing annotation:

- **Blue blocks** mean exon coordinates were present.
- **Orange hatched blocks** mean only a transcript span was available and exon
  boundaries are unknown.
- Orange tracks add the legend
  `Transcript span only (exons unavailable)`.

For example, `gencode.v38.2wayconspseudos.gtf.gz` contains transcript records
but no exon records. Its transcripts can be plotted, but absent exon–intron
structures cannot be reconstructed and are therefore shown as orange hatched
spans.

## 7. Editable PDF output

PDF fonts are embedded as TrueType/font type 42. Adobe Illustrator can edit
labels as text and ungroup or modify the vector components. The plotter does
not embed a raster image.

Example PDFs:

- [Isoform mode](docs/examples/isoform_example.pdf)
- [Reference-region mode](docs/examples/reference_example.pdf)
- Additional generated validation files may be found in `validation_outputs/`.

## 8. Testing and repository layout

Run the automated tests:

```bash
python -m pytest
```

```text
transcript_isoform_plotter/
  cli.py          Command-line interface
  gtf.py          Plain/gzipped GTF parsing
  plotting.py     Editable vector PDF rendering
  splice.py       Splice annotation parsing and matching
tests/            Automated tests
examples/         Example input tables
docs/images/      Optional raster previews
docs/examples/    Editable PDF examples
validation_outputs/ Generated validation PDFs
pyproject.toml    Python package configuration
environment.yml  Conda environment
```

## 9. TODO

- Add color-blind-accessible styles for protein-coding genes, pseudogenes,
  miRNAs, lncRNAs, and other noncoding RNA classes.
- Add a reference-mode legend for gene classes.
- Add exclusion filters for selected biotypes.

## 10. License

MIT
