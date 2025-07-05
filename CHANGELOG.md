# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.1] - 2025-07-05

### Changed
- **Plot Improvements**: Fixed axis units from per m³ to per g solvent for all moment variables
- **Better Readability**: Changed time axis from seconds to minutes across all plots
- **Scientific Notation**: Added proper formatting for large numbers with adequate margins
- **Clean Titles**: Removed technical variable names (y1, y2, etc.) from subplot titles
- **Enhanced Labels**: Improved legend labels to be more descriptive

### Documentation
- **Enhanced README**: Comprehensive installation instructions for VS Code and command line
- **VS Code Integration**: Detailed step-by-step VS Code workflow with Julia extension
- **Installation Verification**: Clear steps to verify successful package installation
- **Better Organization**: Separated installation methods for different user preferences

### Fixed
- Package structure with proper `src/` directory layout
- Precompilation issues - package now installs and imports correctly
- Added result printing for Koopman and Monte Carlo solutions

## [1.0.0] - 2025-07-04

### Added
- Initial stable release of KoopmanCrystallization.jl
- Process modeling with population balance equations
- Deterministic optimization framework
- Koopman operator methods for uncertainty quantification  
- Backoff optimization strategies
- Monte Carlo validation methods
- Complete case study example
- Comprehensive plotting and visualization tools

### Features
- Temperature control optimization for crystallization processes
- Constraint violation tracking and penalty methods
- Root-finding based backoff strategies
- Stochastic analysis using Koopman operators
- Support for multiple uncertainty parameters

### Dependencies
- ModelingToolkit.jl for symbolic modeling
- DifferentialEquations.jl for ODE solving
- Optimization.jl for numerical optimization
- SciMLExpectations.jl for Koopman methods
- Additional scientific computing packages
