# Changelog
All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
### Changed
### Removed
### Fixed

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
