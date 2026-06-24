# Phylogenetic context / outgroup sequences

Drop **per-segment** FASTA files here (or in any directory you point `--context_dir` at) to fold
external reference strains and an outgroup into ARBOR's per-segment trees. The pipeline aligns each
file into **only its own segment** via MAFFT `--add`, so an S-segment reference never pollutes the M/L
trees.

## File names (exact)

| File      | Folded into segment        |
|-----------|----------------------------|
| `S.fasta` | small (`NC_014395_S`)      |
| `M.fasta` | medium (`NC_014396_M`)     |
| `L.fasta` | large (`NC_014397_L`)      |

The suffix is matched to the part after the last `_` in each `--segments` record ID, so this maps
correctly even if you change the reference. A segment with **no** matching file is aligned plain — you
can supply just `S.fasta` and leave M/L empty.

Each FASTA record header becomes a **tip label** in that segment's tree, so name them readably, e.g.
`>ZH548_S`, `>MP12_S`, `>Smithburn_S`.

## Turning it on

Off by default. Enable per run:

```bash
nextflow run tdoerks/arbor --input samplesheet.csv --outdir results \
    --context_dir assets/context
```

## Suggested set for RVFV / MP-12 work

ARBOR maps against the RVFV RefSeq (`NC_014395/6/7`, ZH-548-derived), so the reference already anchors
the centre of each tree. The most informative additions are:

- **MP-12** (the live-attenuated vaccine strain — what these samples *are*): lets you see drift /
  consistency of your preps against the deposited vaccine sequence.
- **A divergent wild RVFV lineage** (e.g. Smithburn, or a Kenya 2006–07 isolate) per segment: gives
  each tree a genuine **outgroup to root on**.

Pull the exact per-segment accessions from **NCBI Virus** (Rift Valley fever virus → segment →
complete sequence) once the strain set is confirmed, concatenate each segment's records into the
matching `S/M/L.fasta`, and re-run with `--context_dir`.
