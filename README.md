# KoopmanCrystallization.jl

A Julia package demonstrating the crystallization process optimization case study presented in the associated research paper (currently under review). It implements Koopman expectation and Backoff-assisted Iterative Optimization strategies for handling parametric uncertainty.


## Installation

Since this is a development package, clone the repository and activate the environment.

### Option 1: Command Line Installation

```bash
git clone https://github.com/Danijel-Zadravec/KoopmanCrystallization.jl
cd KoopmanCrystallization.jl
julia --project=.
```

In Julia REPL:
```julia
using Pkg
Pkg.instantiate()  # Install all dependencies
using KoopmanCrystallization  # Load the package
```

### Option 2: VS Code Workflow (Recommended)

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Danijel-Zadravec/KoopmanCrystallization.jl
   ```

2. **Open in VS Code:**
   - Open VS Code
   - File → Open Folder → Select the `KoopmanCrystallization.jl` folder
   - Install the Julia extension if not already installed

3. **Activate the environment:**
   - Press Ctrl+Shift+P and type "Julia: Start REPL"
   - This will open the Julia REPL in VS Code and should automatically activate the project environment
   - **Verify the environment is active:** Look for `(KoopmanCrystallization) pkg>` in the REPL prompt when you press `]`
   - If not automatically activated, type `] activate .` in the Julia REPL
   - In the Julia REPL, run:
     ```julia
     using Pkg
     Pkg.instantiate()
     ```

4. **Verify installation:**
   ```julia
   using KoopmanCrystallization
   # If no errors, the package is ready to use!
   ```

5. **Run the example:**
   - Open `case_study_final_script.jl` in VS Code
   - Use Ctrl+Shift+P → "Julia: Execute File" or run line by line

## Quick Start

### Running the Case Study

The main example is in `case_study_final_script.jl`. This script demonstrates the complete workflow:

**In VS Code:**
1. Open `case_study_final_script.jl`
2. Make sure Julia REPL is started (Ctrl+Shift+P → "Julia: Start REPL")
3. Run the entire script with Ctrl+Shift+P → "Julia: Execute File"
4. Or execute sections line by line using Ctrl+Enter

**From Command Line:**
```julia
julia --project=. case_study_final_script.jl
```

### What the Example Demonstrates:
- **Model setup and parameterization** - Crystallization process with temperature control
- **Deterministic optimization** - Finding optimal temperature profiles
- **Uncertainty quantification** - Handling parameter uncertainty
- **Backoff optimization** - Constraint violation control using Koopman expectation
- **Results visualization** - Plotting state trajectories and printing key metrics

### Expected Output:
The script will print:
- Koopman Solution Results (objective value and violation frequencies)
- Monte Carlo Results (average metrics from uncertainty analysis)
- Generate plots showing temperature profiles and state evolution


## License

MIT License (see LICENSE file)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use this package in your research, please cite:

