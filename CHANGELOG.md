# Changelog
All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- [Add option to exclude MHC region for CELLECT-MAGMA](https://github.com/perslab/CELLECT/issues/52)
### Changed
- CELLECT-MAGMA gene coordinates updated to extactly match CELLECT-LDSC (which remains the same). See [#51 - Unify LDSC and MAGMA gene coordination files](https://github.com/perslab/CELLECT/issues/51)
### Removed
- Removed `WINDOW_LD_BASED` option (CELLECT-LDSC only feature). 
### Fixed

## [1.2.0] - 2020-04-16
### Added
- CELLECT-MAGMA (analysis types supported: prioritization, conditional). See https://ctg.cncr.nl/software/magma for more information about the algorithm and auxiliary files. 
- Extra output directory layer. The results are output into 2 separate subdirectories of the base output directory: CELLECT-LDSC and CELLECT-MAGMA respectively.
- WiKi (respective sections for CELLECT-MAGMA)
### Changed
- Minimized the overlap of functions between CELLECT-LDSC and CELLECT-MAGMA via includes. 
- Config file. The file was divided into sections of common and LDSC-/MAGMA-specific parameters. 
- WiKi (respective sections for CELLECT-LDSC)
- README Documentation
## [1.1.0] - 2020-04-14
### Added
### Changed
- Improved CELLECT-LDSC handling of continuous gene annotations. Importantly, this changes the prioritization results slightly. 
CELLECT v. 1.1.0 and 1.0.0 produce prioritization results with Pearson's correlation ~0.95.
### Removed
- Redundant "make multigeneset" functionality
### Fixed
- Some legacy references to incorrect paths
### Technical details
See commit for overview of important changes: https://github.com/perslab/CELLECT/commit/c7513a19c0dfc02b905e507f4c94924cee6a4ca4
See also https://github.com/perslab/CELLECT/issues/8 for a description of the problem.
- Replaced use of BEDtools for merging overlapping genes into a single track with BEDOPs which finds all unique overlapping regions and is more appropriate for continuous data
- General restructuring of rules and scripts to account for above
- Added BEDOPs to env
- Some additional tidying and minor bugfixes
## [1.0.0] - 2020-01-01
### Added
- CELLECT-LDSC conditional analysis type
- CELLECT-LDSC h2 analysis type
- Matrix expression specificity input file format
- Result parser (results/*.csv files)
- README documentation
- log files in Snakemake workflow
### Changed
- Result file headers
- Config file restructured (breaking changes)
### Removed
- Several config settings
### Fixed
- Snakemake workflow optimization (wildcards etc)


## [0.1.1] - 2019-08-20
First stable release for CELLECT S-LDSC prioritization. No support for S-LDSC h2 or conditional analysis.
