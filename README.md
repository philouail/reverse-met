# Reverse Metabolomics: a generic reanalysis workflow for public LC-MS/MS data

A reproducible end-to-end **workflow for public-data reanalysis at scale**,
built on the [RforMassSpectrometry](https://www.rformassspectrometry.org/)
ecosystem. It takes one or more target compounds, finds them across public LC-MS/MS
deposits by spectral library search, confirms each hit with SIRIUS, and then
measures the confirmed compounds at MS1 in every biological file of every
confirmed assay (for MassIVE, every file a metadata source describes plus every
file with a confirmed hit), including the many files where the instrument never
fragmented them.

The workflow is generic in the targets and in their number. Two applications
ship with it, under `application/`, each with a pair:

| application | pair | question |
|---|---|---|
| `omeprazole` | omeprazole / 5-hydroxyomeprazole (exogenous) | CYP2C19 metaboliser activity: a genotyped, dosed cohort as ground truth, then the ratio in blood across the other deposits |
| `kyn-trp` | L-tryptophan / L-kynurenine (endogenous) | the kynurenine/tryptophan ratio against age in adult plasma and serum, across cohorts |

![The workflow: six numbered steps, the packages of the thesis each step uses (hexes) and the external tools it calls (grey)](docs/workflow.png)

*The workflow. Steps 1 to 5 are the same for every application; step 6 is where an
application asks its own question. Hexes mark the packages developed in the thesis,
and grey names the external tools a step calls.*

Each step is one notebook in `application/<name>/`:

| step | notebook | what it does |
|---|---|---|
| 1 Search | `01-search.qmd` | resolves each target to its library spectra (FASSTrecords, by the first InChIKey block) and their matches (cosine of at least 0.7 over at least 3 peaks), sets aside library spectra whose precursor m/z is not a plausible ion of the target (no adduct within 25 ppm), and keeps the deposits that carry every target |
| 2 Metadata | `02-metadata.qmd` | first removes the deposits whose ReDU record already rules them out (another species; in kyn-trp also no blood); then retrieves each deposit's sample metadata from MassIVE (Pan-ReDU, plus the submitter's own table, which it looks for among the deposit's files), MetaboLights (ISA-Tab) and the Metabolomics Workbench (mwTab), maps it onto common fields, and assigns each hit to its assay |
| 3 Curation | `03-curation.qmd` | removes deposits that cannot answer the question, each with a recorded reason, keeps the assays with a hit for every target, and removes a MassIVE deposit whose files mix set-ups |
| 4 Confirmation | `04-confirmation.qmd` (+ `confirmation/*.qmd`) | re-scores every hit with SIRIUS and CSI:FingerID and grades it; keeps the assays in which every target is confirmed, less those set aside by hand with a reason; sets, per assay and target, the ion to extract, its m/z and the range of confirmed retention times, from the best grade the target reached; reads each assay's instrument from the headers of its raw files, which sets the m/z tolerance of step 5 |
| 5 Extraction | `05-extraction.qmd` | measures each target at MS1 in every biological file of each confirmed assay: ion trace over step 4's range and a margin beyond it, at the m/z tolerance of the assay's instrument, a retention-time window set on the MS1 peaks the confirmed spectra point to, peak picking, fixed peak criteria and the isobar rule |
| 6 Analysis | `06-analysis.qmd` | the application's question, one notebook per application: one unit per person and matrix, one assay kept per set of samples, then the application's own models |

What differs between the applications is their `compound.csv`, their hand
corrections and checks in `curation/`,
omeprazole's hand-obtained `dataset_metadata/` files, the parts of steps 1 to 4
marked **"Manual work"**, the text about each application's own targets and
deposits, and the code that goes with them, omeprazole's extra notebook `_05-probe-checks.qmd` (a child of
`05-extraction.qmd`), and step 6. `05-extraction.qmd` and the three confirmation
notebooks are identical between the applications;
`06-analysis.qmd` is written for each application, and the two share their helpers
in `R/analysis.R`.

---

## What this is

This repo serves two intertwined purposes:

**1. A methods showcase.** The workflow shows how
[`RuSirius`](https://github.com/RforMassSpectrometry/RuSirius),
[`Chromatograms`](https://www.bioconductor.org/packages/Chromatograms/),
and the Bioconductor MS backends
([`MsBackendMassIVE`](https://www.bioconductor.org/packages/MsBackendMassIVE/),
[`MsBackendMetaboLights`](https://www.bioconductor.org/packages/MsBackendMetaboLights/),
[`MsBackendMetabolomicsWorkbench`](https://www.bioconductor.org/packages/MsBackendMetabolomicsWorkbench/))
compose into one pipeline that:

- harmonises sample-level metadata across three repository schemas
  (Pan-ReDU, ISA-Tab, mwTab) with an auditable trail (step 2);
- runs SIRIUS-based formula and structure confirmation at scale via
  `RuSirius`, per repository (step 4);
- extracts targeted ion traces with `Chromatograms` over `Spectra` objects,
  measures each peak with `MsQuality`, and produces one row per file with each
  target's detection and peak area across every assay (step 5).

**2. Two questions.** Omeprazole is a proton-pump inhibitor cleared primarily by
hepatic CYP2C19, and the ratio of parent to 5-hydroxyomeprazole is a long-standing
proxy for CYP2C19 phenotype, whose distribution varies strikingly by ancestry. The
`omeprazole` application first tests the pipeline on a dosed, genotyped cohort,
where the answer is known, then asks whether the ratio in blood moves with
anything the other deposits record. The `kyn-trp` application asks whether the
kynurenine/tryptophan ratio, a read-out of indoleamine 2,3-dioxygenase activity,
rises with age in adult blood across cohorts: an endogenous pair where presence
says little and the ratio carries the question.

### Why MS1 extraction, and not MS2 alone

This is the pipeline's central methodological claim.

These deposits are acquired by **data-dependent acquisition (DDA)**: the
instrument decides *during the run* which ions to fragment, choosing
mostly by intensity, within a limited duty cycle. A compound can be
present at a perfectly measurable level and receive **no MS2 scan at
all**, simply because the instrument was busy elsewhere. So an MS2-only
search (MASST and relatives) cannot distinguish *"the compound is
absent"* from *"the compound is present but was never fragmented"*, and
because selection is intensity-driven, the misses are biased toward the
**low-abundance** cases. **MS2-only prevalence is a biased lower bound.**

The two levels are complementary: **MS2 gives identity** (fragments
constrain the structure, and can separate isomers); **MS1 gives
coverage** (recorded continuously, for everything, with no selection
step). The workflow uses MS2 *once* to pin the mass and retention time,
then measures at MS1 in **every** file.

**This has a direct precedent in proteomics: *match-between-runs*.**
MBR identifies a peptide by MS2 in the run where it happened to be
fragmented, then transfers that identification to other runs by accurate
mass and aligned retention time, quantifying at MS1 where the precursor
was never selected. **Metabolomics repository reanalysis has largely not
adopted the equivalent step**, and that gap is what this pipeline fills.

Two consequences worth stating plainly:

- **It makes cross-repository comparison fairer, not just larger.**
  Comparing deposits by MS2 detection rate substantially compares
  *instruments and acquisition settings* (a deeper precursor list
  fragments more, so it "detects" more). MS1 is acquired the same way
  regardless, so extraction removes that confound.
- **It trades specificity for sensitivity, and that cost is real.** MS2
  can separate some compounds that share a mass, but not one that gives the
  same product ions; MS1 cannot separate any. Extraction extends coverage
  only as far as chromatography preserves identity. Omeprazole sulfone (the
  CYP3A4 product) is an exact isomer of 5-hydroxyomeprazole, and
  fragmentation does not tell them apart. It elutes near omeprazole. Step 5
  keeps it out in three ways. On reverse phase the peak choice is weighted
  toward the centre of the window. The window is kept modest. Because the
  `isobar` column of `compound.csv` names the sulfone, `resolve_isobars()`
  also drops a metabolite peak lying nearer omeprazole's retention time than
  its own, each retention time being the median of the MS1 apexes the
  confirmed spectra point to. The retention-time weighting is the main
  defence: without it, a window wide enough to reach the sulfone measures the
  wrong compound. Omeprazole has an isomer of its own, a hydroxylated
  omeprazole sulfide, which `compound.csv` names too. Where most MS2 spectra
  at omeprazole's mass are of that isomer, step 4 sets the assay aside by hand
  (MSV000088255).

---

## Current state

> [!NOTE]
> The pipeline runs in the six steps above, with one step-6
> notebook per application. Every number is computed in the notebooks. Each step
> also writes its summary to `artifacts/summary_step<N>.csv`.

The audit trail for every curation decision is in the notebooks and in `R/`,
each marked as manual work or as a rule by hand: the deposits step 2 removes by hand, the constants,
column maps and related deposits in `02-metadata.qmd`, the Workbench file order in
`R/metadata.R` (`MWB_FILE_ORDER`), the exclusion list and the stray hits in
`03-curation.qmd`, the assays set aside in `04-confirmation.qmd` (`SET_ASIDE`),
the peak corrections in `curation/`, and the rule for kyn-trp's sleep study in its
`06-analysis.qmd`. The inputs those decisions act on are mostly regenerated rather
than versioned: the metadata cached under `application/<name>/dataset_metadata/`
and the FASSTrecords extracts under `application/<name>/StructureMASST_output/`
(see below). There are two exceptions in `dataset_metadata/`. Six files that
omeprazole's step 2 obtained by hand are tracked. The Workbench files
`ST*_AN*.txt` are neither tracked nor fetched by any code, yet step 5 reads them
(step 5 skips an assay
without its file). See [dataset_metadata/](#dataset_metadata-per-application)
below.

---

## Setup

1. **R and packages.** R 4.6 with Bioconductor 3.23, Quarto 1.5 or later, and the
   packages listed under [Dependencies](#dependencies). Every notebook checks the
   pinned versions when it starts (`R/setup.R`).
2. **SIRIUS** 6.3 or later on the `PATH`: step 4 calls it through `RuSirius`. Log in
   once through the SIRIUS GUI; the notebooks use the login it caches.
3. **FASSTrecords.** Download `masst_records.zip` (42.6 GB) from Zenodo,
   [doi:10.5281/zenodo.18199544](https://doi.org/10.5281/zenodo.18199544), and
   unzip it on a disk with room for the database it holds, `masst_records.sqlite`
   (about 141 GB). Tell R where it is in `~/.Renviron`:

   ```
   FASST_DB=/path/to/masst_records.sqlite
   ```

   Steps 1 and 2 query it and write their extracts to
   `application/<name>/StructureMASST_output/`.
4. **Disk space for the raw files.** The SIRIUS runs, step 4's checks and step 5
   download the raw files of the assays they read through the `MsBackend`
   packages, into one BiocFileCache.

## Running it

1. **Set up** R, SIRIUS and the FASSTrecords database (see [Setup](#setup)).
2. **Pick an application** and work from its directory: `application/omeprazole/`
   or `application/kyn-trp/`. Nothing renders from the repo root.
3. **Render the steps in order**, `01` to `06`, one at a time.

```r
setwd("application/kyn-trp")   # or application/omeprazole

quarto::quarto_render("01-search.qmd")
quarto::quarto_render("02-metadata.qmd")
quarto::quarto_render("03-curation.qmd")
# Step 4 combines the outputs of the three confirmation notebooks
# (artifacts/hits-confirmed-*.csv, git-ignored). A plain render stops if none
# exist and reuses them as they are if they do. On a fresh clone, or after steps
# 1-3 change, run SIRIUS first (one instance at a time; the three repository
# notebooks run one after the other):
#   quarto render 04-confirmation.qmd -P run_confirmations:true
quarto::quarto_render("04-confirmation.qmd")
quarto::quarto_render("05-extraction.qmd")
quarto::quarto_render("06-analysis.qmd")
```

Or from the repo root, with retries on transient repository failures:

```bash
bash bench/rerun.sh kyn-trp 01 02 03 massive metabolights workbench 04 05 06
```

Leave out `massive metabolights workbench` to reuse existing confirmation
outputs. With no stage named, it renders `04 05 06`.

The SIRIUS runs, step 5 and omeprazole's step 4 (its MS2 check of MSV000088255)
download the raw files through the relevant `MsBackend` (large; allow time and
disk), all into one BiocFileCache. That cache takes one writer at a time: two
renders at once deadlock, so render one notebook at a time. `bench/rerun.sh` runs
its stages one after the other for this reason. SIRIUS authenticates with the
login cached by its GUI; log in there first.

The first render of `02` fills `dataset_metadata/`: the Pan-ReDU tables
(`<id>_panredu.tsv`), each MassIVE deposit's own sample table, found among the
deposit's table files (`<id>_metadata.tsv`), the MetaboLights investigation and
assay files (`<id>_i_Investigation.txt`, `a_*.txt`), the study files (`s_<id>.txt`), and each Workbench deposit's
metadata and file listing as `MsBackendMetabolomicsWorkbench` returns them
(`<id>_mwb_metadata.rds`, `<id>_mwb_files.rds`; delete one to fetch it again).
Step 2's first pass reads ReDU from the
FASSTrecords database and keeps the rows it uses in
`StructureMASST_output/redu_files.csv`. No code fetches the Workbench analysis
files `<id>_AN*.txt`; the ones step 5 reads are kept in git.

Step 5 skips any MetaboLights or Workbench assay whose `s_` or `_AN` file is
missing.

`R/setup.R`, sourced at the top of every notebook, stops if one of five pinned packages is missing or
off-version (RuSirius exactly 1.0.4; MsBackendMassIVE >= 0.99.1,
MsBackendMetaboLights >= 1.7.4, MsBackendMetabolomicsWorkbench >= 0.99.0,
Chromatograms >= 1.3.3), or if the `sirius` binary is not on PATH or is older
than 6.3. It does not check any other package. See the note on
MsBackendMetaboLights under [Dependencies](#dependencies).

---

## Appendix A · Repository layout

### Notebooks

Every notebook lives under `application/<name>/`; paths below are relative to
that directory. Every step notebook also writes its summary,
`artifacts/summary_step<N>.csv`.

| notebook | step | reads | writes |
|---|---|---|---|
| `01-search.qmd` | 1 | `compound.csv`, the FASSTrecords database (path in `FASST_DB`, default `E:/FASSTrecords/db/masst_records.sqlite`) or its cached extract `StructureMASST_output/fasstrecords_hits.csv` | `artifacts/hits_search.csv`; also, when the FASSTrecords database is present and `StructureMASST_output/fasstrecords_hits.csv` is absent or was queried for other targets, that file with `fasstrecords_targets.csv` (the InChIKeys it was queried for) and `fasstrecords_deposits.csv` (each target's deposits with a hit) |
| `02-metadata.qmd` | 2 | `artifacts/hits_search.csv`, the FASSTrecords database or `StructureMASST_output/redu_files.csv`, `dataset_metadata/*` | `artifacts/hits_assigned.csv`, `artifacts/assay_metadata.csv`, `artifacts/metadata_entities.csv`, `artifacts/related_deposits.csv`, `artifacts/removed_step2.csv` (+ `artifacts/metadata_inventory.csv` in kyn-trp, `artifacts/rescue_samples.csv` in omeprazole); `StructureMASST_output/redu_files.csv` when the database is present; and the `dataset_metadata/` cache: Pan-ReDU `*_panredu.tsv`, MassIVE submitter tables `MSV*_metadata.tsv`, MetaboLights `*_i_Investigation.txt`, `a_*.txt` and `s_*.txt`, Workbench `ST*_mwb_metadata.rds` and `ST*_mwb_files.rds`, and in omeprazole `MSV000082493_metadata.tsv` and `MSV000084008_metadata.tsv` |
| `03-curation.qmd` | 3 | `artifacts/hits_assigned.csv`, `artifacts/metadata_entities.csv`, `artifacts/assay_metadata.csv`, `artifacts/removed_step2.csv` (+ `artifacts/metadata_inventory.csv` in kyn-trp) | `artifacts/hits_curated.csv`, `artifacts/dataset_metadata_table.csv`, `artifacts/exclusions.csv` (+ `figures/metadata_inventory_kept.{png,pdf}` in kyn-trp) |
| `confirmation/massIVE-hit-confirmation.qmd` | 4 | `artifacts/hits_curated.csv` | `artifacts/hits-confirmed-massive.csv`, `results/MSV*/`, `sirius_projects/MSV*/` |
| `confirmation/metaboLights-confirmation.qmd` | 4 | `artifacts/hits_curated.csv` | `artifacts/hits-confirmed-metaboLights.csv`, `results/MTBLS*/`, `sirius_projects/MTBLS*/` |
| `confirmation/metabolomics-workbench-confirmation.qmd` | 4 | `artifacts/hits_curated.csv` | `artifacts/hits-confirmed-metabolomics-workbench.csv`, `results/ST*_AN*/`, `sirius_projects/ST*_AN*/` (one per assay) |
| `04-confirmation.qmd` | 4 | `artifacts/hits-confirmed-*.csv`, `artifacts/dataset_metadata_table.csv`, `artifacts/hits_curated.csv`, `artifacts/assay_metadata.csv`, `results/<id>/<role>_results.csv` (one per target), and the headers of raw files with a confirmed spectrum (already on disk from the confirmation notebooks, found by name in the BiocFileCache the MsBackend packages download into; nothing is downloaded for them). Omeprazole also reads every MS2 scan at omeprazole's mass in MSV000088255's sample files (`load_spectra_ms2()`, which downloads any file not yet in the cache), once: the classes are kept in `artifacts/msv88255_ms2_check.csv` | `artifacts/hits-confirmed-all.csv`, `artifacts/datasets-confirmed.csv`, `artifacts/extraction_windows.csv` (per assay and target: ion, m/z and range of confirmed retention times), `artifacts/file_instruments.csv` (per file whose header was read: manufacturer, model, mass analyser, ion source and polarity; cached, so a render reads only files not yet in it), `artifacts/assay_instruments.csv` (per assay: instrument, class orbitrap / tof / unknown, the m/z tolerance step 5 uses, polarity of its files and how many differ, and what the metadata records beside it) |
| `05-extraction.qmd` | 5 | `artifacts/hits-confirmed-all.csv`, `artifacts/datasets-confirmed.csv`, `artifacts/extraction_windows.csv`, `artifacts/assay_instruments.csv` (each assay's m/z tolerance), `artifacts/dataset_metadata_table.csv`, `dataset_metadata/*` (MetaboLights `s_<deposit>.txt` and Workbench `<ST>_AN*.txt` are required, and an assay without its file is skipped; the MassIVE tables are read through `msv_native_sample_table()`), `curation/*` (optional; `curation/tolerance_check.csv` lists the assays checked again), `artifacts/hits_curated.csv` (only for an assay under `tolerance_check`), `peaks_review_v2/_extract_times.csv` (optional, runtime table only) | `artifacts/bio_{files,assays,meta}_per_ds.rds`, `artifacts/per_dataset/<dataset_id>.rds`, `artifacts/ms1_windows.csv` (the window in which each target's peak was looked for, per assay), `artifacts/tolerance_check_<dataset_id>.rds` (per assay under `tolerance_check`), `artifacts/sample_table_perfile.csv`, `artifacts/sample_table_all.csv`, `artifacts/positive_control.csv`, `figures/step5_windows.{png,pdf}`, `figures/step5_uncertain_windows.{png,pdf}` (where a window is uncertain), `figures/step5_pc_misses.{png,pdf}` |
| `_05-probe-checks.qmd` (omeprazole only, on the dosed cohort: the offset of MS2 trigger times from the MS1 apex, the sulfide fragment check of the parent channel, each extraction choice turned off in turn with the genotype model refitted and the pre-dose negatives counted, and the window-pad sweep toward the sulfone) | 5 | child of `05`, so it also uses 05's in-memory `results`, `conf_final`, `bio_files_per_ds` and `ms1_windows`. Files: `curation/rt_apex_overrides.csv`, `curation/ms1_drop_list.csv`, `artifacts/hits_curated.csv`, `dataset_metadata/MSV000084008_metadata.tsv`, `dataset_metadata/MSV000082493_metadata.tsv`, `peaks_review_v2/_cache/MSV000084008.rds`, and the raw MS2 of the confirmed parent files (only when `artifacts/isobar_fragment_check.csv` is missing) | `artifacts/isobar_fragment_check.csv` (cached) |
| `06-analysis.qmd` (kyn-trp) | 6 | `artifacts/sample_table_perfile.csv`, `artifacts/bio_files_per_ds.rds`, `artifacts/bio_meta_per_ds.rds`, `artifacts/datasets-confirmed.csv`, `artifacts/hits-confirmed-all.csv`, `artifacts/related_deposits.csv`, `artifacts/metadata_entities.csv`, and the submitter tables `dataset_metadata/<deposit>_metadata.tsv` of MSV000082630 (the sleep study's labels), MSV000096564, MSV000084556 and MSV000094097 (the first-time-point check) | `artifacts/assay_choice.csv`, `artifacts/ratio_covariate_pooled.csv`, `artifacts/ratio_covariate_per_assay.csv`, `figures/ratio_vs_age_by_matrix.png`, `figures/ratio_vs_age_healthy.{png,pdf}`, `figures/loo_assay_by_matrix.{png,pdf}` |
| `06-analysis.qmd` (omeprazole) | 6 | `artifacts/sample_table_perfile.csv`, `artifacts/bio_files_per_ds.rds`, `artifacts/bio_meta_per_ds.rds`, `artifacts/datasets-confirmed.csv`, `artifacts/hits-confirmed-all.csv`, `artifacts/related_deposits.csv`, and the dosed study's tables `dataset_metadata/MSV000084008_metadata.tsv` and `dataset_metadata/MSV000082493_metadata.tsv` | `artifacts/assay_choice.csv`, `artifacts/ratio_scan.csv`, `figures/ratio_volunteer.png` |

### Applications (`application/<name>/`)

The two applications, `application/omeprazole/` and `application/kyn-trp/`, each
run the workflow for one pair and are self-contained (each has its own
`compound.csv`, `dataset_metadata/`, `artifacts/`,
`confirmation/`, `curation/`, `StructureMASST_output/`, `results/`,
`sirius_projects/`, `peaks_review_v2/`). The repo root holds only what the
applications share, with the licence and citation files: the `R/` library,
the `bench/` tooling and the workflow figure (`docs/workflow.png`).

| file | what it is |
|---|---|
| `compound.csv` | the targets, one row per target: role, name, formula, [M+H]+ m/z and InChIKey first block, and the optional `isobar` and `isobar_at` columns (see *The targets* below). Both applications list a `parent` and a `metabolite` |

**The targets.** `compound.csv` has one row per target, any number from one, and
its row order is the target order of every loop, table and column. A role is a
unique lowercase label, a letter then letters or digits, and it names files and
directories (`results/<id>/<role>_results.csv`, `sirius_projects/<id>/<role>/`).
`load_compound_spec()` refuses the names in `RESERVED_ROLES`: `par` and `met`,
names already given to fields or labels (`all`, `none`, `both`, `neither`,
`metab`, `total`, `ppm`) and the Windows device names (`con`, `nul`, `com1`,
...). A target's columns and fields carry a prefix (`target_prefix()`): `par`
for the parent, `met` for the metabolite, the role itself for any other
target, so a third target `serotonin` adds `serotonin_auc` beside `par_auc`
and `met_auc`. Every step asks for every target: step 1 keeps the deposits with
a hit of each, step 3 the assays with a hit of each, step 4 the assays where
each is confirmed, and step 5 counts a file as confirmed by MS2, or detected in
MS1, only when each target is.

`isobar` names a compound that shares a target's m/z, and `isobar_at` the target
at whose retention time it elutes (`isobar_reference()`). With two targets
`isobar_at` may be left blank, as in both applications, and is then the other
target; otherwise loading `compound.csv` stops on an isobar without one.

The step summaries word the targets by their number: "deposits with both
targets" with two, "deposits with the target" with one, "deposits with all 3
targets" with three. The rejected assays listed as "... one target confirmed"
with two targets are "... some but not all targets confirmed" with any other
number. The detection tiers depend on the number too (*Detection tiers* in
Appendix B). A few diagnostics compare a pair and run only with two targets:
step 5's *Checking the tolerance* and omeprazole's probe checks, which need the parent and the metabolite. Step 6
is written for its application's pair.

Application figures land in `application/<name>/figures/`.

### R helper library

Reusable functions; each notebook sources only the files it needs.

| file | contents |
|---|---|
| `R/setup.R` | the version check of the pinned packages and the SIRIUS binary (see [Dependencies](#dependencies)), then `library()` calls for the shared core packages (Spectra, the three MsBackends, xcms, RuSirius, dplyr/tidyr/ggplot2, ...); creates `artifacts/`; the helpers that name the targets: `target_prefix()` / `role_of_prefix()` (a role's column prefix and back) and `targets_phrase()` / `some_targets_phrase()` (the step summaries' wording). Step 1 attaches RSQLite itself when it builds the hit extract. Step 6 does not source this file and attaches dplyr and ggplot2 itself, and the confirmation notebooks load their own packages instead of sourcing it |
| `R/masst.R` | `infer_adduct()` and `ADDUCT_CANDIDATES` (step 1's library-spectrum check); `load_compound_spec()`, which checks the application's `compound.csv` (`RESERVED_ROLES`) and builds `compound_spec` + `role_formula` from it; `isobar_reference()` (the target each declared isobar elutes at); `results_dir_for()` and `structure_ranks()` (where a target ranks among the structures of its formula); `PHENO`, the CYP2C19 activity order omeprazole's validation regresses on |
| `R/metadata.R` | metadata retrieval and harmonisation: `redu_first_pass()` (step 2's check of what ReDU already records), `fetch_panredu_metadata()`, `fetch_massive_submitter_table()` (finds a MassIVE deposit's own sample table among its table files and caches it), `msv_native_sample_table()` (Pan-ReDU plus submitter columns joined by filename stem; a submitter-only file is admitted if `massive_list_files()` lists it and it carries at least `min_fields` canonical fields, 2 by default and 0 where step 5 reads the sample type), `native_sample_table()`, `harvest_metadata_coverage()`, `fetch_assay_metadata()`, `load_metadata()`, `assay_lookup_for()` / `file_to_assay()` (which assay a hit's file belongs to), `MWB_FILE_ORDER` / `mwb_analysis_order()` / `mwb_listed_files()` (the Workbench file-to-analysis link, declared by hand where the default order is wrong), `select_bio_files_massive()`, `EXCL_ST` / `EXCL_NONBIO` / `is_nonbio_file()` (blanks, QCs and standards, by declared sample type and by file name), `norm_file_key()` / `to_bare_name()`, `is_missing_val()`, `short_ds()`, `mwb_cached()` / `massive_files()` (the Workbench and MassIVE listings, read once and kept in `dataset_metadata/`), `is_copy()` (the files a deposit holds again from a related deposit, counted once, from step 3 on) |
| `R/ms1.R` | steps 4 and 5. The peak picker and its settings: `walk_bound()`, `pick()`, `gate_reason()` / `accept_peak()`, `add_msquality()`, `snap_one()` and `acceptance_window()` (the window), `is_dropped()`. Spectra: `load_spectra()`, `load_spectra_ms1()` (step 5). The instrument: `PPM_BY_CLASS` / `PPM_DEFAULT` / `instrument_class()`, `cached_raw_paths()` and `read_instrument()` (step 4 reads the headers). Step 4's ion, m/z and range: `dominant_adduct()`, `best_tier()`, `derive_extraction_config()`, `window_setters()`, and `extraction_window()` (step 5 reads them back). The traces: `extract_eics()`, and for the tolerance check `scan_centroids()` and `common_ions()`. Per assay: `per_target()` (one field of every target, by prefix), `process_eics()`, `resolve_isobars()`, `build_sample_table()`, `curation_inputs()` / `curation_key()`, `anchor_key()`, `PICKER_VERSION`, `repick_cached()`, `dataset_meta_row()` / `row_is_hilic()`, `run_ms1_for_dataset()`. The ablation helpers that re-run the extraction from the cached traces with one choice changed are in omeprazole's `_05-probe-checks.qmd`. Peak helpers: `count_peaks()`, `flag_fwhm()`, `check_confirmed()` |
| `R/analysis.R` | step 6, shared by the two `06-analysis.qmd`: `has_metadata()` and `ms2_confirmed()` (counted as step 5 counts), `same_sample_groups()` / `keep_one_assay()` (one assay per set of samples), `unitise()` / `disagree_text()` (one unit per person and matrix), `canon_site()`, `harmonise_sex()`, `nz()` |
| `R/peak_review_panels.R` | manual-audit tool, run per application (`Rscript R/peak_review_panels.R <application> [dataset_id ...]`): re-picks from the cached wide traces and writes per-assay review panels (each rejection labelled with the criterion behind it, also in `index.csv`) and the extraction-timing log to `application/<name>/peaks_review_v2/` |

### bench/

`bench/rerun.sh` renders the steps of one application in order, retrying when a
repository does not answer. Every analysis lives in the notebooks.

| script | what it does |
|---|---|
| `bench/rerun.sh` | re-renders one application's stages in the order named, retrying transient repository failures (`bash bench/rerun.sh <application> [stage ...]`, stages `01` to `06`, and `massive` / `metabolights` / `workbench` for one confirmation notebook; `04 05 06` when no stage is named). A stage the application has no notebook for is skipped. Stages run one after the other, since the BiocFileCache is single-writer and two renders at once deadlock. Confirmation is uncapped unless `CONFIRM_CAP` is set. Logs to `bench/rerun-<application>-<notebook>.log` |

### artifacts/ (per application)

Each application writes its intermediate files to `application/<name>/artifacts/`.
The whole tree is **derived and git-ignored**: every file is regenerated by
re-running the notebook that writes it (see the notebook table above), in step
order. Nothing here is a primary input.

### curation/ (per application)

Hand-made corrections to step 5, made after looking at the peak panels; the
step-5 notebook marks them as manual work. All four files are optional: without one,
an application extracts without that kind of correction.

| file | what it records |
|---|---|
| `curation/ms1_drop_list.csv` | peaks to drop, keyed by (`dataset`, `role`, `file`), the file matched by `norm_file_key()`, each with a stated reason. A dropped peak's gate reads `manual drop` |
| `curation/ms1_accept_list.csv` | peaks to accept, keyed by (`dataset`, `role`, `file`) and the `apex_rt` read off the peak panel, each with a stated reason: the file is picked within 5 s of that apex instead of in the window, and the peak criteria still apply. An accepted peak's gate reads `manual accept`; the positive control counts these apart (`by_hand`) |
| `curation/tolerance_check.csv` | the assays step 5 extracts once more at the time-of-flight tolerance, to see what a wider tolerance would add (*Checking the tolerance*); without the file the check is skipped |
| `curation/rt_apex_overrides.csv` | the MS1 apex retention time to anchor the window on where a confirmed MS2 scan's automatic snap landed on the wrong peak. Rows are matched by (dataset_id, role) and ms2_rt within 2 s, each with a note; the override replaces that scan's anchor, and the picker still chooses each file's peak inside the resulting window |

kyn-trp drops both targets' peaks in one file of MSV000085607 (a distorted
injection the deposit ran again) and has no overrides file. Omeprazole's drops are
all in ST002044_AN003327: peaks of the shifted compound that step 5's *Checking
the tolerance* found next to each target. Its overrides file holds only its
header. A deposit's corrections are part of its cache key (`curation_key()`), so
adding or removing one re-picks that assay from its stored traces. Read by
`R/ms1.R` (extraction), `R/peak_review_panels.R` (audit panels) and omeprazole's
`_05-probe-checks.qmd`.

### dataset_metadata/ (per application)

Cached native metadata files, **git-ignored**. Step 2 and `R/metadata.R`
re-fetch the Pan-ReDU tables (from the live API, not a pinned version), the
MetaboLights investigation, study and assay files, and each MassIVE deposit's own sample
table. Omeprazole's step 2 also downloads `MSV000082493_metadata.tsv` from its
authors' repository at a pinned commit and derives `MSV000084008_metadata.tsv`
from it. Six files obtained by hand for omeprazole are kept in git. One other
group is neither fetched by any code nor kept in git, yet step 5 needs it: the
Workbench mwTab files (`ST*_AN*.txt`), without which step 5 skips the assay.

| pattern | source | how it gets cached |
|---|---|---|
| `MSV*_metadata.{tsv,csv}` | submitter sample tables | found by `fetch_massive_submitter_table()` in step 2: any table file of the deposit, in any folder, with a column that names the deposit's own MS files, downloaded through `MsBackendMassIVE`. Only the rows of the deposit's files are kept, and the table is cached as `MSV*_metadata.tsv`. A deposit searched without result leaves `MSV*_metadata_all.none` and is not searched again; nothing is cached when a download fails. A table already present is used as it is: omeprazole's step 2 writes MSV000082493 (from the authors' GitHub repository at commit 95586620) and MSV000084008 (its plasma rows), and MSV000086975 was obtained by hand and is tracked |
| `MSV*_panredu.tsv` | GNPS2 Pan-ReDU | fetched by `fetch_panredu_metadata()` from `https://redu.gnps2.org/...`; an empty `*_panredu.none` marks an accession Pan-ReDU has no rows for |
| `MSV*_panredu.json` | GNPS2 Pan-ReDU, exported by hand (omeprazole only: MSV000085256, MSV000082433) | no code can fetch it again; tracked in git (not ignored) and parsed by omeprazole's `02-metadata.qmd` (`pj()`) |
| `MTBLS*_i_Investigation.txt` | MetaboLights ISA-Tab investigation file | downloaded from the EBI FTP by `mtbls_inv_assay_names()` |
| `a_MTBLS*_*.txt` | MetaboLights ISA-Tab assay file | cached by `fetch_assay_table()` from `mtbls_assay_data()` (`MsBackendMetaboLights`), falling back to the EBI FTP if that fails; the two `a_MTBLS12231_*.txt` files are tracked |
| `s_MTBLS*.txt` | MetaboLights ISA-Tab study file | cached by `native_sample_table()` (step 2) from `mtbls_sample_data()` on first read; step 5 (`05-extraction.qmd`) skips an MTBLS assay whose `s_` file is missing. `s_MTBLS12231.txt` is tracked |
| `ST*_AN*.txt` | Workbench mwTab analysis file | not fetched by the pipeline: placed in `dataset_metadata/` by hand (git-ignored) and read only by step 5 (`05-extraction.qmd`), which skips a Workbench assay whose file is missing. Step 2 reads the same metadata through `MsBackendMetabolomicsWorkbench::mwb_metadata()` and caches it as `<id>_mwb_metadata.rds` |

**Two-source policy for MassIVE deposits.** The Pan-ReDU table is the base, and
the submitter's extra columns are added to its rows by filename stem. A
submitter row that Pan-ReDU lacks is added only if the deposit's inventory
(`massive_list_files()`) lists the file and at least two canonical fields carry
a value. Step 5 reads the declared sample type on every row, however thin
(`msv_native_sample_table(min_fields = 0)`), since the sample type alone says a
file is a blank, QC or standard. Which files are extracted is decided from the
inventory in step 5 (`select_bio_files_massive()`). ReDU/MIxS missing-value
placeholders (`"missing value"`, `"not applicable"`, ...) are filtered everywhere
via `is_missing_val()`. The rationale and audit findings live in **Appendix C**
below.

### results/ and sirius_projects/ (per application)

Both are produced locally by the confirmation notebooks and are **git-ignored**.
Re-running them needs SIRIUS 6.3 or later and the CSI:FingerID web service.

```
sirius_projects/<id>/<role>/batchNNN/<role>_bNNN.sirius   # batched (default, CONFIRM_BATCH=100): one SIRIUS project per batch
sirius_projects/<id>/<role>/batchNNN_results.csv         # that batch's checkpointed results
sirius_projects/<id>/<role>/<role>.sirius                # CONFIRM_BATCH=Inf: one SIRIUS project per role
results/<id>/{hits_rt,<role>_results}.csv               # SIRIUS output, read by step 4
```

Here `<role>` is a target's role in `compound.csv`: `parent` or `metabolite` in
both applications.

`<id>` is the deposit for MassIVE and MetaboLights and the assay (deposit and
analysis, e.g. `ST002044_AN003325`) for the Workbench, whose deposits often hold
several assays of the same samples in separate files.

The confirmation notebooks take three environment switches: `CONFIRM_ONLY`
(comma-separated deposits to process, MassIVE only), `CONFIRM_CAP` (MassIVE and
Workbench only, off by default: send only the top-N hits by cosine per role and
biofluid in MassIVE and per role in the Workbench. It is a quick anchor check,
not the confirmation the analysis uses. The MetaboLights notebook is always
uncapped.) and `CONFIRM_BATCH` (annotate each role in batches of N features, 100
by default; `Inf` runs a role in one piece). A batch is its own SIRIUS
project with its results checkpointed beside it, at
`sirius_projects/<id>/<role>/batchNNN_results.csv`, and a later render reuses
every finished batch, so an outage or a restart costs only the batch in flight.
Batching changes durability, not results: no step couples features, so a
feature's annotation does not depend on which batch it lands in. A render stops
if stored results no longer match the retention-time table they are numbered
against, rather than attach them to the wrong spectra.

### StructureMASST_output/ (per application)

The FASSTrecords extracts, all regenerated from the database and none of them in git:

- `fasstrecords_hits.csv`, the hit set: written by the FASSTrecords chunk of
  `01-search.qmd` and read by step 1 otherwise. The chunk runs only when the
  database is present and the file is absent or was queried for other targets
  than `compound.csv` lists. As step 1 writes it,
  it holds only the hits in the deposits that carry every target. Without the
  database, a change of targets stops step 1.
- `fasstrecords_targets.csv` (omeprazole only), written with the hit set: the InChIKeys it was
  queried for, which step 1 compares with `compound.csv`, so renaming a target
  needs no new query. A hit set written before this file is checked on the
  InChIKeys of its hits.
- `fasstrecords_deposits.csv`, written with the hit set: each target's deposits
  with a hit over the whole archive (`inchikey_first_block`, `deposit`), with the
  same check of the library spectra as the hits. Step 1's "deposits with a hit"
  counts them. Without this file it counts the deposits in the hit set, which,
  for a hit set limited to the deposits with every target, are only those.
- `redu_files.csv`: ReDU's per-file rows (species, body site, sample type) for
  the deposits with a hit, written by step 2's first pass from the database's
  ReDU table when the database is present, and read by it otherwise.

A fresh clone therefore needs either the FASSTrecords database (DOI
10.5281/zenodo.18199544, path in `FASST_DB`), from which steps 1 and 2 regenerate
these files, or `fasstrecords_hits.csv` and `redu_files.csv` supplied separately
and placed here; without them, steps 1 and 2 fail.

---

## Appendix B · Key concepts

**Confidence grades** (SIRIUS, per confirmed feature; rules from each confirmation notebook):

| grade | rule |
|---|---|
| **strict** | formula rank 1 AND structure rank 1 |
| **standard** | formula rank ≤ 3 AND structure rank ≤ 5 |
| **liberal** | formula rank ≤ 15 AND structure rank ≤ 15 |
| **not_confirmed** | beyond liberal; dropped before `hits-confirmed-*.csv` |

All three confirmed grades flow through step 4 and are preserved in
`hits-confirmed-all.csv`; the `confidence` column lets the analysis stratify or
weight by grade. Step 4 drops a confirmation from a blank, QC or standard file
(`is_nonbio_file()`). It sets each assay's ion, m/z and range from the
confirmations at the best grade the target reached there (`best_tier()`), while
step 5 anchors its window on confirmations of every grade.

**What step 4 sets for step 5** (`04-confirmation.qmd`, `R/ms1.R`):

- per assay and target (`derive_extraction_config()`, in
  `extraction_windows.csv`): the ion, the molecular ion if it was confirmed and
  otherwise the most frequent confirmed adduct; its m/z, the median observed
  precursor m/z of that ion; and the range, from the earliest to the latest
  confirmed retention time of that ion, padded by `RT_PAD_S` = 5 s each side;
- per assay, the instrument, read from the headers of its raw files with a
  confirmed spectrum (`read_instrument()`): five of those files, spread over their
  sorted list (`N_HEADER`), or all of them where the files read disagree on instrument
  class or polarity. Only files already on disk are read (`cached_raw_paths()`),
  and each file read is kept in `file_instruments.csv`. The assay takes the class
  its files record (`instrument_class()`), `unknown` where they record none or
  more than one, and the m/z tolerance of that class (`assay_instruments.csv`);
- assays set aside by hand: `SET_ASIDE` in `04-confirmation.qmd` lists, under
  manual work, assays that confirmed every target by the rules but are set aside
  with a reason; they leave everything step 4 writes. Omeprazole sets aside
  MSV000088255, where most MS2 spectra at omeprazole's mass are of its isomer, a
  hydroxylated omeprazole sulfide (a chunk classes every MS2 scan at omeprazole's
  [M+H]+ by marker fragments). kyn-trp sets aside none.

**How step 5 measures a target** (`R/ms1.R`; `05-extraction.qmd` explains each
choice):

| setting | value |
|---|---|
| m/z tolerance | per assay, from step 4: 5 ppm on an Orbitrap, 20 ppm on a time-of-flight, 20 ppm where the headers do not tell (`PPM_BY_CLASS`, `PPM_DEFAULT`) |
| trace read | step 4's range, plus `WIDE_EXTRA` = 50 s each side; it bounds the trace, not the window |
| anchor | each confirmed MS2 scan moved to the apex of the MS1 peak it sits on, in its own file (`snap_one()`: the flanks are walked as the picker walks them), and dropped if more than `SNAP_W` = 20 s away; a hand override replaces it |
| one anchor per file | the file's most intense anchored apex, so a file fragmented many times weighs as one |
| window | the best-supported cluster of anchors (`CLUST_TOL` = 10 s) padded by `APEX_PAD` = 5 s each side, then extended over accepted peaks within `COHERE_S` = 1 s of its edge, up to 2 × `APEX_PAD` from the cluster. Where two clusters are equally supported it spans every anchor, both clusters included, and is flagged uncertain; with no anchor it falls back to step 4's range |
| peak choice | reverse phase: area weighted by a Gaussian on the apex's distance from the window centre (`RT_SIGMA` = 4 s); HILIC: largest area |
| peak criteria | an apex of at least 500 counts, SNR of at least 12 (HILIC 3), FWHM of 1 to 30 s (HILIC at least 1.3 s); HILIC also needs a return to baseline (prominence at least 0.25 of the apex) and no neighbour more than 1.3 times taller within 15 s; a 3-point reverse-phase peak needs that prominence too |
| truncated peak | a trace that empties before half height, shaped like one peak (its rise and fall may dip by up to `TRUNC_DIP` = 20% of the apex): it can be detected, but has no area |
| isobar rule | only for a target whose `compound.csv` row names an isobar (`resolve_isobars()`): a detected peak at least as near the retention time of the target in `isobar_at` (with two targets, the other one by default) as its own is rejected, each retention time being the median of the anchors its window rests on. It stands down where either of the two windows is uncertain or their centres lie less than `CHAN_SEP` = 6 s apart |
| `auc0` | each peak's area also on a baseline started from 0 at a boundary beside an empty scan; recorded for step 6's checks, never used to detect or choose a peak |

The picker is at version 5 (`PICKER_VERSION`). Each assay's extraction is cached
in `artifacts/per_dataset/<dataset_id>.rds`. A change to step 4's ion, m/z or
range, or to the assay's m/z tolerance, reads the raw files again. A change of
picker version, of chromatography gate, of the confirmed anchors (`anchor_key()`)
or of the assay's hand corrections (`curation_key()`) re-picks the assay from its
stored traces, reading no raw file (`repick_cached()`).

Step 5 measures every biological file of an assay. Blanks, QCs and standards are
removed by the sample type the metadata declares (`EXCL_ST`), read on every row of
a MassIVE table, and by one narrow file-name rule applied to every repository
(`EXCL_NONBIO`: a name starting with "pool", or with "std" and a digit, or a "qc" or "blank"
token). A confirmed file whose MS1 peak is missed is reported in the positive
control, with the rule behind each miss.

*Checking the tolerance*, which runs only with two targets, extracts each assay
listed in `curation/tolerance_check.csv` again at 20 ppm, with the same
ion, m/z, range, picker and hand corrections, and reads its raw centroids.
`scan_centroids()` gives, per MS1 scan of a target's window, the most intense
centroid within 20 ppm: its deviation and intensity. A file is *shifted* when the wider tolerance adds it to the files with
both targets and both its apexes lie beyond the assay's own tolerance.
`common_ions()` finds the ions present in at least 90% of the files over the
parent's window and 15 s either side (centroids of at least 100,000 counts, in
0.01 m/z bins), and each file's calibration offset is the median deviation of
those ions from their median m/z in the files that are not shifted. The result is
cached in `artifacts/tolerance_check_<dataset_id>.rds`.

**How step 6 counts** (`R/analysis.R` and each `06-analysis.qmd`):

- only files a metadata source covers are analysed (`has_metadata()`);
- one assay per set of samples (`keep_one_assay()`): the assays of one deposit
  that share at least half the sample names of the smaller, and the panels step 2
  declares measured on two set-ups (`same_samples_two_assays`). The kept assay is
  the one where MS1 detects both targets in the most files, then the one with
  more files where MS2 confirms both, then the first assay id, counted over every
  extracted file. A study deposited twice (`same_study_different_arm`) is one set
  of files, not a set of assays;
- one unit per person and matrix (`unitise()`), or one per file where no subject
  is named. A unit's ratio is the median over its files that give one. Its age,
  sex, health and other fields take the value all its files agree on, and are
  missing where the files disagree; those disagreements are reported
  (`disagree_text()`);
- in kyn-trp, health in three levels (`health_level()`): healthy for a value in a list of
  healthy words, never a pattern; not recorded for a missing value or a
  placeholder; not healthy otherwise.

kyn-trp fits the log2 kynurenine/tryptophan ratio against age in adult (18 or
older) blood units with a sex recorded, one model per matrix (plasma, serum), with a term for each
assay and for sex and health where they vary inside an assay. A model is fitted
where it has at least 8 units spanning at least 10 years. Each matrix is fitted
again on its healthy units alone, and each assay alone. One rule is written by
hand: in the sleep study MSV000082630, only the files its metadata labels healthy
are kept. Two checks stay in the notebook, outside the chapter's results: the
models on the first time point of each person, ordered only from the deposit's
own metadata, never from file names; and the models on the `auc0` areas.

Omeprazole tests the dosed cohort (MSV000084008 plasma, MSV000082493 skin). The
negative controls are the plasma draws at 0 and 5 min and the forehead skin
swabbed before the first dose, each with an exact binomial upper bound on the
rate of false detections. Each volunteer's ratio is that of the two targets'
areas under the curve over the draw times (trapezoids, a target not detected
counted as zero), per study day, then averaged over days. It is regressed on the
CYP2C19 metaboliser class (`PHENO`, one step per class), and refitted without each
volunteer in turn. The rest of the corpus is scanned: every metadata field of the
per-file table, tested against the ratio in blood inside each assay where it
varies (F test), all p values corrected together by Benjamini-Hochberg at a 5%
false discovery rate. Four checks stay in the notebook, outside the chapter's
results: the first study day alone, a single draw near 180 min, sex, age and
matrix tested across assays, and the model on the `auc0` areas.

**Detection tiers** (per file, set in step 5 by `build_sample_table()`). With two
targets, as in both applications, the first row of `compound.csv` is the ratio's
numerator:

| tier | first target's peak | second target's peak | ratio (`log2_ratio`) |
|---|---|---|---|
| `both` | yes | yes | log2 of the first target's peak area over the second's, `log2(par/met)` here; none where a peak is truncated |
| `<role>_only` (`parent_only`) | yes | no | none (right-censored) |
| `<role>_only` (`metab_only`, shortened for the metabolite) | no | yes | none (left-censored) |
| `neither` | no | no | none; the file contributes nothing |

With one target, or three or more, the tier is `all`, `none` or the roles detected
joined by `+` (`parent+serotonin`), and there is no `log2_ratio` column.

Step 6 builds a ratio only where both targets are detected with an area, and
imputes nothing. The one exception is the dosed cohort's areas under the curve,
where a target not detected in a draw counts as zero.

**Exclusion vocabulary** (`artifacts/removed_step2.csv` from step 2, and
`artifacts/exclusions.csv` from step 3, which carries step 2's too):

| reason | step | meaning |
|---|---|---|
| `not_human` | 2 | no file ReDU describes is human, and at least one is of another species (judged only where ReDU describes every file with a hit); kyn-trp also lists one deposit by hand |
| `no_blood` | 2 (kyn-trp) | every file ReDU describes is from a body site other than blood |
| `not_reachable` | 2 | a repository no MsBackend package reads |
| `metadata_unreadable` | 2 (kyn-trp) | the deposit's list of assays could not be read |
| `duplicate_deposit` | 3 (kyn-trp, by hand) | the same files deposited again under another accession |
| `no_or_limited_metadata` | 3 (omeprazole, by hand) | no usable biological metadata anywhere (deposit, Pan-ReDU, publication), or only a single value for every field |
| `access_restricted` | 3 (omeprazole, by hand) | demographics behind a controlled-access portal |
| `dia` | 3 (omeprazole, by hand) | data-independent acquisition; single-scan MASST matching is unreliable |
| `in_vitro` | 3 (omeprazole, by hand) | cell culture or reference-library deposit |
| `non_human` | 3 (omeprazole, by hand) | non-human samples, where ReDU carries no taxonomy to catch them |
| `reference_material` | 3 (omeprazole, by hand) | a reference-material suite, not a study cohort |
| `metadata_unreachable`, `no_metadata`, `not_human`, `no_blood_found`, `no_adult_with_age` | 3 (kyn-trp) | from step 2's inventory: metadata that could not be read when step 2 ran, none found, another species, no blood, or no adult with an age |
| `no_assay_with_both_roles`, `no_assay_with_every_target` | 3 | no assay of the deposit has a hit of every target: their hits land in different assays (the second code with any other number of targets than two) |
| `no_hit_matched_to_assay` | 3 | step 2 matched none of the deposit's hits to an assay |
| `mixed_setups` | 3 | a MassIVE deposit whose files mix set-ups (chromatography or polarity) |

Assays dropped at SIRIUS confirmation (step 4) are not in this list: they appear
in `hits_curated.csv`, and step 4 lists them with the formula SIRIUS ranked first
for each role, and the assays it sets aside by hand with their reason.

---

## Appendix C · Metadata sources for MassIVE deposits

A MassIVE deposit can be described by two sources: the table its submitter
uploaded, and [Pan-ReDU](https://redu.gnps2.org/), which re-harmonises submitter
metadata against ReDU/MIxS vocabularies. Neither is complete, and they disagree in
ways that decide which files are measured and which biology can be recovered. What
we found while building steps 2 and 5:

- **Pan-ReDU can miss files.** On MSV000082493 it did not index the blood
  collection added in a later deposit update, together with that deposit's
  CYP2C19 columns, the ones the omeprazole validation depends on. Omeprazole's
  step 2 reads that deposit's metadata from its authors' repository instead
  (manual work).
- **Submitter tables vary in completeness, naming and place.** Some advertise
  `.mzXML` where MassIVE hosts `.mzML` (MSV000088255); one lists files that do
  not exist on the server at all (MSV000095143, excluded in step 3). A table can
  sit in any folder of the deposit, or in a later update, so step 2 looks at every
  table file and keeps one that names the deposit's own MS files. The Chagas
  serum panel (MSV000092754, MSV000092756, MSV000092757 and MSV000092758) names its subject and serology columns in
  its own way, so step 2 maps them by hand.
- **Controlled-vocabulary fields often hold placeholders** (`"missing value"`,
  `"not applicable"`), filtered everywhere by `is_missing_val()`.

The policy the code implements:

1. **Join.** Pan-ReDU is the base; the submitter's extra columns are added to its
   rows by filename stem. A submitter row Pan-ReDU lacks is admitted only if the
   deposit's own inventory (`massive_list_files()`) lists the file and at least two
   canonical fields carry a value (`msv_native_sample_table()`).
2. **Which files are measured** is decided in step 5 from the deposit's inventory,
   not from either metadata source: every file a metadata source describes, plus
   every file with a confirmed hit (`select_bio_files_massive()`), unless the
   metadata declares that file a blank, QC or standard. Blanks, QCs and standards
   are removed by the declared sample-type column, read on every row, and a narrow
   file-name rule.
3. **Two accessions sharing file names** are either one study deposited twice (the
   same files: count once) or one panel measured twice (different data: keep
   apart). Only reading the deposits tells the two apart, so step 2 declares each
   case by hand in `related_deposits.csv`, and step 6 keeps one assay of a panel
   measured twice.

That a public-data pipeline must reconcile three independent, individually
incomplete things (spectral hits, the file inventory and every metadata source)
is itself a finding, separate from the MS1-recovery claim.

---

## Dependencies

```r
BiocManager::install(c(
  "Spectra",
  "xcms", "MsCoreUtils", "MsQuality", "MetaboCoreUtils",
  "Chromatograms", "mzR", "BiocFileCache"
))
install.packages(c("dplyr", "tidyr", "ggplot2", "jsonlite", "curl", "readxl",
                   "RSQLite", "patchwork"))
remotes::install_github("RforMassSpectrometry/RuSirius@v1.0.4")
# The gabri branches keep the assay in the cached file name, so a deposit's
# positive- and negative-mode files with the same name do not overwrite each other.
remotes::install_github("RforMassSpectrometry/MsBackendMassIVE@gabri")
remotes::install_github("RforMassSpectrometry/MsBackendMetaboLights@gabri")
remotes::install_github("RforMassSpectrometry/MsBackendMetabolomicsWorkbench")
```

`R/setup.R` (sourced at the top of every notebook) enforces:

| package | constraint |
|---|---|
| RuSirius | exactly 1.0.4 (JOSS submission pin) |
| MsBackendMassIVE | ≥ 0.99.1 |
| MsBackendMetaboLights | ≥ 1.7.4 |
| MsBackendMetabolomicsWorkbench | ≥ 0.99.0 |
| Chromatograms | ≥ 1.3.3 |
| SIRIUS binary | ≥ 6.3 |

**MsBackendMetaboLights.** The pipeline expects the fix for assay IDs (1.7.9
and later) and the merged MetaboLights metadata (1.7.10 and later), while
`R/setup.R` accepts 1.7.4. That version is safe
for these data, because only the positive-mode assays of MTBLS11656 and MTBLS1866
were downloaded. Installing the fix renames the cache files of MTBLS11656, so the
next extraction that reads its raw files downloads them once more.

---

## Licence and citation

The code in this repository is released under the **Apache License 2.0**
(see [`LICENSE`](LICENSE)). The public mass-spectrometry data it reanalyses
remains under the terms set by each original depositor; the accessions are
listed in `artifacts/datasets-confirmed.csv` and in the chapter.

If you use this workflow or its results, cite it with the metadata in
[`CITATION.cff`](CITATION.cff).
