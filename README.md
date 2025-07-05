# KoopmanCrystallization.jl

A Julia package demonstrating the crystallization process optimization case study presented in the associated research paper Optimal Control of Uncertain Batch Processes via Koopman Expectation-Assisted Gradient Methods (currently under review). It implements Koopman expectation and Backoff-assisted Iterative Optimization strategies for handling parametric uncertainty.


## Installation & Quick Start

### Prerequisites
- Julia 1.6 or higher
- VS Code with Julia extension (recommended)

### Option 1: VS Code Workflow (Recommended)

1. **Clone and open:**
   ```bash
   git clone https://github.com/Danijel-Zadravec/KoopmanCrystallization.jl
   ```
   - File → Open Folder → Select `KoopmanCrystallization.jl`

2. **Setup environment:**
   - Ctrl+Shift+P → "Julia: Start REPL" (auto-activates project environment)
   - In Julia REPL:
     ```julia
     using Pkg
     Pkg.instantiate()
     using KoopmanCrystallization  # Verify installation
     ```

3. **Run the example:**
   - Open `case_study_final_script.jl`
   - Ctrl+Shift+P → "Julia: Execute File" or execute line by line with Shift+Enter

### Option 2: Command Line Workflow

1. **Clone and setup:**
   ```bash
   git clone https://github.com/Danijel-Zadravec/KoopmanCrystallization.jl
   cd KoopmanCrystallization.jl
   julia --project=.
   ```

2. **Install dependencies:**
   ```julia
   using Pkg
   Pkg.instantiate()
   using KoopmanCrystallization  # Verify installation
   ```

3. **Run the example:**
   ```bash
   julia --project=. case_study_final_script.jl
   ```

### Expected Output
- Koopman Solution Results (objective and violation frequencies)
- Monte Carlo Results (uncertainty analysis averages)
- Plots of temperature profiles and state evolution


## License

MIT License (see LICENSE file)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use this package in your research, please cite:

