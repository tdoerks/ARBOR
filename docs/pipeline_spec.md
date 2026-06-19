# arbor

**ARBOR** — **A**mplicon **R**econstruction, **B**ootstrapping & **O**utbreak **R**esolution. A viral
**tiling-amplicon** pipeline: Illumina paired-end reads → primer-trimmed mappings → low-frequency
variants (iSNVs) + per-segment consensus → core phylogeny. A Nextflow/nf-core port of the *CLC Genomics
Workbench* RVFV (Rift Valley Fever Virus, MP-12, tri-segment S/M/L) SOP.

Reusable for any tiling-amplicon virus (give it a reference + primer BED); RVFV multi-segment handling
is the default because the reference is multi-record.

GitHub repo: **ARBOR** · internal nf-core name: `arbor`

## Overview

```
reads ─► FastQC ─► fastp ─► Bowtie2(--local) ─► sort+index ─┬─► QC (samtools stats, mosdepth)
                                                            ▼
                                          iVar trim (primer BED) ─► sort+index ─┬─► QC (post-trim)
                                                                                ▼
                          ┌──────────────────── whole-reference BAM ───────────────────┐
                          ▼                                                             ▼
              iVar variants  +  LoFreq        (iSNVs)            split BAM by segment (S/M/L)
                          │                                                  ▼
                          ▼                                      iVar consensus (per segment)
                       VCF/TSV                                               ▼
                                                     group by segment ─► MAFFT(--add context) ─► IQ-TREE
                                                                                                    ▼
                                                                                       per-segment trees
```

Variant calling runs on the whole multi-segment reference. **Consensus + phylogeny are per-segment**
(segmented viruses reassort — S/M/L can have different trees), so the BAM is split by segment before
consensus, regrouped across samples, aligned (with optional external reference strains), and a tree is
built per segment.

## Samplesheet schema

CSV header `sample,fastq_1,fastq_2`. Paired-end only.

| Column | Required | Notes |
|---|---|---|
| `sample` | yes | unique id → `meta.id` |
| `fastq_1` | yes | R1 `.fastq.gz` |
| `fastq_2` | yes | R2 `.fastq.gz` |

`meta.single_end = false` (the tiling/primer-trim path assumes PE).

## Steps

| # | Process | Module | I/O (in → out) | nf-core module |
|---|---|---|---|---|
| 1 | Read QC | `FASTQC` | `[meta, reads]` → zip/html | `nf-core/fastqc` |
| 2 | Adapter/quality trim | `FASTP` | `[meta, reads, adapter?]` → `reads`, `json` | `nf-core/fastp` |
| 3 | Reference index | `BOWTIE2_BUILD` | `[meta, ref_fasta]` → `index` | `nf-core/bowtie2/build` |
| 3b | Reference faidx | `SAMTOOLS_FAIDX` | `[meta, ref_fasta]` → `fai` | `nf-core/samtools/faidx` |
| 4 | Map (local) | `BOWTIE2_ALIGN` | `[meta, reads] + index + ref + save_unaligned + sort_bam` → sorted `bam` | `nf-core/bowtie2/align` |
| 5 | Index BAM | `SAMTOOLS_INDEX` | `[meta, bam]` → `bai` | `nf-core/samtools/index` |
| 6 | Mapping QC | `SAMTOOLS_STATS`, `MOSDEPTH` | `[meta, bam, bai] (+ref/bed)` → stats, coverage | `nf-core/samtools/stats`, `nf-core/mosdepth` |
| 7 | **Primer trim** | `IVAR_TRIM` | `[meta, bam, bai] + primer_bed` → trimmed `bam` | `nf-core/ivar/trim` |
| 8 | Re-sort + index | `SAMTOOLS_SORT`, `SAMTOOLS_INDEX` | `[meta, bam]` → sorted `bam` + `bai` | `nf-core/samtools/sort`, `…/index` |
| 9 | Post-trim QC | `SAMTOOLS_STATS`, `MOSDEPTH` | (as #6) | (same) |
| 10 | Variants (freq) | `IVAR_VARIANTS` | `[meta, bam] + ref + fai + gff + save_mpileup` → `tsv` | `nf-core/ivar/variants` |
| 11 | Variants (low-freq) | `LOFREQ_CALL` (+`LOFREQ_FILTER`) | `[meta, bam, bai] + ref + fai` → `vcf` | `nf-core/lofreq/call`, `…/filter` |
| 12 | Split BAM by segment | `SAMTOOLS_VIEW` | `[meta, bam, bai] + ref`, region=segment (ext.args) → per-segment `bam` | `nf-core/samtools/view` |
| 13 | Consensus (per segment) | `IVAR_CONSENSUS` | `[meta+seg, bam] + ref + save_mpileup` → `fasta` | `nf-core/ivar/consensus` |
| 14 | Merge consensus | `CAT_CAT` | `[seg, [fastas…]]` → one fasta/segment | `nf-core/cat/cat` |
| 15 | Select segment recs | `SEQKIT_GREP` | `[meta, seqs] + pattern` → per-segment fasta | `nf-core/seqkit/grep` |
| 16 | MSA (+context) | `MAFFT_ALIGN` | `[seg, fasta] + [ctx, context_fasta]` (`--add`) → alignment | `nf-core/mafft/align` |
| 17 | Phylogeny (per segment) | `IQTREE` | `[seg, aln, []]` → `treefile` | `nf-core/iqtree` |
| 18 | Report | `MULTIQC` | fastqc + fastp + samtools + mosdepth + ivar → report | `nf-core/multiqc` |

**Every tool has a free nf-core module — zero custom modules.**

## Data flow (channel wiring)

Shared reference handled as **value channels** (one reference for the whole run):
```
ch_ref     = [ [id:'ref'], file(params.reference) ]            // value
BOWTIE2_BUILD(ch_ref) -> ch_index   (value, .first())
SAMTOOLS_FAIDX(ch_ref) -> ch_fai    (value, .first())
ch_primer_bed = file(params.primer_bed)                        // bare path value
ch_context    = params.context_fasta ? [[id:'ctx'], file(params.context_fasta)] : [[],[]]
```

Per-sample mapping & trim:
```
FASTQC(reads); FASTP(reads, adapter) -> ch_trim
BOWTIE2_ALIGN(ch_trim, ch_index, ch_ref, false, true) -> ch_bam   // sort_bam=true
SAMTOOLS_INDEX(ch_bam) -> bai ; join -> [meta,bam,bai]
// REVIEW: bowtie2 emits .csi not .bai; ivar/trim needs bai -> SAMTOOLS_INDEX explicitly
SAMTOOLS_STATS(...); MOSDEPTH(... amplicon/primer bed ...)               // pre-trim QC

IVAR_TRIM([meta,bam,bai], ch_primer_bed) -> trimmed bam (unsorted)
SAMTOOLS_SORT -> SAMTOOLS_INDEX -> ch_final = [meta, bam, bai]
SAMTOOLS_STATS(ch_final); MOSDEPTH(ch_final)                            // post-trim QC
```

Variants (whole reference):
```
IVAR_VARIANTS(ch_final.map{m,b,i->[m,b]}, ch_ref_fa, ch_fai, [], false) -> tsv
LOFREQ_CALL(ch_final, ch_ref_fa, ch_fai) -> vcf ; LOFREQ_FILTER -> filtered vcf
// REVIEW: lofreq may need indelqual/alnqual preprocessing — confirm at build
```

Per-segment consensus + phylogeny (the fan-out/fan-in — the riskiest wiring):
```
ch_seg = channel.fromList(params.segments)        // e.g. ['S','M','L'] = reference record IDs
ch_seg_bam = ch_final.combine(ch_seg)             // [meta,bam,bai,seg]
SAMTOOLS_VIEW(region=seg via ext.args) -> [meta+[seg], seg_bam]
IVAR_CONSENSUS([meta+seg, seg_bam], ch_ref_fa, false) -> [meta+seg, consensus.fa]
  .map { meta, fa -> [ meta.seg, fa ] }.groupTuple()       // group across samples by segment
  -> CAT_CAT -> [ [id:seg], combined.fa ]
MAFFT_ALIGN([id:seg, combined.fa], ch_context /*--add*/, [[],[]],[[],[]],[[],[]],[[],[]], false) -> aln
IQTREE([id:seg, aln, []], [],[]…) -> per-segment treefile
// REVIEW: per-segment fan-out (combine) then fan-in (groupTuple by seg) — verify on real RVFV data
// REVIEW: segment IDs in params.segments MUST match the reference FASTA record names exactly
```

MultiQC: `.mix()` FASTQC.zip, FASTP.json, SAMTOOLS_STATS.stats (pre+post), MOSDEPTH summaries,
IVAR_VARIANTS/TRIM logs (all natively MultiQC-supported).

## Parameters

| Param | Default | Configurable? | Source |
|---|---|---|---|
| `input` | — | required | samplesheet |
| `outdir` | — | required | — |
| `reference` | — | required path | RVFV S/M/L FASTA |
| `primer_bed` | — | required path | tiling primer BED |
| `segments` | `['S','M','L']` | configurable list | reference record IDs |
| `context_fasta` | null | configurable path | external NCBI strains for tree |
| `adapter_fasta` | null (auto-detect) | configurable | fastp; CLC uses Illumina universal adapters |
| `fastp_qualified_quality` | 13 | configurable | CLC quality limit 0.05 ≈ Q13 |
| `fastp_length_required` | 50 | configurable | CLC min length 50 |
| `fastp_n_base_limit` | 2 | configurable | CLC max ambiguities 2 |
| `bowtie2_args` | `--local` | configurable | CLC local alignment, len/sim 0.9 |
| `ivar_trim_args` | `-e -q 20 -m 50` | configurable (`-e` = keep no-primer reads, hardcoded intent) | CLC "uncheck remove whole read" |
| `ivar_variants_min_freq` | 0.01 | configurable (`-t`) | CLC min frequency 1% |
| `ivar_variants_min_depth` | 250 | configurable (`-m`) | CLC min coverage ~250x |
| `ivar_variants_min_qual` | 20 | configurable (`-q`) | CLC central Q≥20 |
| `lofreq_sig` | 0.005 | configurable | CLC significance 0.5% |
| `lofreq_min_cov` | 100 | configurable | CLC min coverage 100–500x |
| `consensus_min_depth` | 10 | configurable (ivar `-m`) | CLC low-coverage threshold 10 |
| `consensus_min_freq` | 0 | configurable (ivar `-t`) | CLC "Vote" majority |
| `consensus_n` | `N` | hardcoded (ivar `-n N`) | CLC undetermined = N |
| `mafft_args` | `--maxiterate 1000 --localpair` | configurable | CLC "Very Accurate" = L-INS-i |
| `iqtree_model` | `MFP` | configurable | CLC Model Testing → ML |
| `iqtree_bootstrap` | 1000 | configurable (`-B`) | CLC bootstrap 100/1000 |
| `slurm_account` | null | configurable | Beocat/Ceres profiles |
| `skip_fastqc` / `skip_lofreq` / `skip_ivar_variants` / `skip_phylogeny` | false | configurable | nf-core convention |

## Resource estimates

Viral genomes are tiny (S/M/L ≈ 12 kb total), so most steps are light. Beocat `resourceLimits` caps apply.

| Process | Label | Rationale |
|---|---|---|
| FastQC, fastp, samtools/*, mosdepth, ivar/* | `process_low` | small viral BAMs, fast |
| BOWTIE2_BUILD / faidx | `process_single` | trivial on a ~12 kb reference |
| BOWTIE2_ALIGN | `process_low` | amplicon depth, small reference |
| LOFREQ_CALL | `process_medium` | per-position model, the heaviest caller here |
| MAFFT_ALIGN | `process_low` | few short sequences per segment |
| IQTREE | `process_medium` | 1000 UFBoot + ModelFinder |

## Profiles

- **`beocat`** — KSU Beocat (reuse the phinder/pangea profile pattern).
- **`ceres_ars`** — USDA Ceres overlay.
- **`test`** — tiny samplesheet + a small reference + primer BED for stub/CI.

## Open questions / decisions locked

- [x] Per-segment phylogeny (S/M/L separate trees) — segmented virus / reassortment. User-confirmed (recommended).
- [x] Both variant callers (iVar variants + LoFreq) output for cross-check. User-confirmed.
- [x] Optional `--context_fasta` external strains folded into the MSA via MAFFT `--add`. User-confirmed.
- [x] Bowtie2 (`--local`) mapper. User-confirmed.
- [ ] `params.segments` record IDs must match the reference FASTA headers exactly → mark `// REVIEW`.
- [ ] LoFreq pre-processing (`indelqual`/`alnqual`) — confirm necessity at build.
- [ ] Optional GFF for iVar variants AA annotation — left empty (`[]`) unless an RVFV GFF is provided.

## Readiness checklist

- [x] Tool list agreed (CLC → iVar amplicon toolchain).
- [x] Data flow validated against module I/O signatures (incl. per-segment fan-out/fan-in).
- [x] Every step maps to an existing nf-core module (zero custom).
- [x] Samplesheet schema defined (PE).
- [x] Configurable vs hardcoded params decided, defaults proposed from the SOP.
- [x] Resource labels assigned.
