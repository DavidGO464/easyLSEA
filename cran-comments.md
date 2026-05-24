## R CMD check results

0 errors | 0 warnings | 1 note (local and GitHub Actions)
0 errors | 0 warnings | 3 notes (win-builder R-devel)

## Test environments

* Local: macOS Tahoe 26.2, R 4.5.1 (aarch64-apple-darwin20)
* GitHub Actions (R-CMD-check.yaml):
  - ubuntu-latest (R release)
  - ubuntu-latest (R devel)
  - ubuntu-latest (R oldrel-1)
  - macOS-latest (R release)
  - windows-latest (R release)
  All 5 platforms: 0 errors | 0 warnings | 1 note
* win-builder (R-devel): 0 errors | 0 warnings | 3 notes (see below)

## Notes

### Note 1 — Non-standard top-level file
'cran-comments.md' is intentionally included for CRAN submission notes.
This is standard practice for CRAN submissions.

### Note 2 — Possibly misspelled words
The following terms appear in DESCRIPTION and may be flagged as possibly
misspelled. They are established terms in the lipidomics and bioinformatics
fields:

* LSEA — Lipid Set Enrichment Analysis; widely used in lipidomics literature
* fgsea — the official name of the Bioconductor package by Korotkevich et al.
  (doi:10.1101/060012); intentionally lowercase to match the package name
* lipidomics — the large-scale study of cellular lipids; standard term in
  the field (cf. lipidomics.net, Journal of Lipid Research)

### Note 3 — Title case for 'fgsea'
win-builder suggests capitalizing 'fgsea' as 'Fgsea' in the Title field.
We intentionally keep the lowercase 'fgsea' to match the official Bioconductor
package name and its published citation. Capitalizing it would be misleading
to users familiar with the package.

### Note 4 — CITATION.cff
'CITATION.cff' is a GitHub citation file (https://citation-file-format.github.io)
added to the repository for GitHub's "Cite this repository" feature.
It has been added to .Rbuildignore and is not included in the package tarball.

## Downstream dependencies

None.
