# KoopmanCrystallization.jl

A Julia package demonstrating the crystallization process optimization case study presented in the associated research paper (currently under review). It implements Koopman expectation and Backoff-assisted Iterative Optimization strategies for handling parametric uncertainty.


## Installation

Since this is a development package, clone the repository and activate the environment:

```bash
git clone https://github.com/yourusername/KoopmanCrystallization.jl.git
cd KoopmanCrystallization.jl
julia --project=.
```

In Julia:
```julia
using Pkg
Pkg.instantiate()
```

## Quick Start


See `case_study_final_script.jl` for a complete example demonstrating:
- Model setup and parameterization
- Deterministic optimization
- Uncertainty quantification
- Backoff optimization with constraint violation control


## License

MIT License (see LICENSE file)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use this package in your research, please cite:

