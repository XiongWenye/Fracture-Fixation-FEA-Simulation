# Fracture Fixation FEA Simulation

A comprehensive finite element analysis simulation toolkit for investigating the biomechanical behavior of bone-fixator systems during fracture healing. This project models stress transfer mechanisms and analyzes how different fixator stiffnesses influence the healing process through computational simulation.

## 🎯 Project Overview

This project utilizes 2D finite element analysis to study the critical biomechanical phenomenon of "stress transfer" during fracture healing. The simulation models:

- **Bone-fixator system** with time-dependent callus maturation
- **Stress redistribution** among bone, fixator, and callus tissues
- **Fracture gap strain evolution** under different fixator rigidities
- **Parametric studies** comparing flexible, standard, and rigid fixators

## 🏗️ Project Structure

```txt
tree
├── main.py                      # Main simulation runner
├── stress_and_model.py          # Stress analysis utilities (original output)
├── mesh.ipynb                   # Jupyter notebook for mesh visualization
├── utils/                       # Core simulation modules
│   ├── fea_core.py              # Finite element analysis core
│   ├── materials.py             # Material property definitions
│   ├── simulation.py            # Simulation workflow
│   ├── analysis.py              # Stress and strain calculations
│   └── plot_utils.py            # Visualization utilities
├── output/                      # Basic simulation results
├── output_advanced/             # Parametric study results
├── assets/                      # Generated animations and visualizations
└── report/                      # LaTeX report and documentation
```

## 🚀 Quick Start

### Prerequisites

```bash
pip install numpy matplotlib scipy pandas jupyter
```

### Running Simulation

```bash
python main.py
```

This will:
1. Run parametric studies for three fixator types (Flexible, Standard, Rigid)
2. Generate stress distribution visualizations
3. Create comparison plots and summary data
4. Save all results to `output_advanced/` directory

## 📊 Key Features

### 1. Material Modeling
- **Bone**: Elastic modulus 18 GPa, Poisson's ratio 0.3
- **Callus**: Time-dependent stiffening from 0.1 GPa to 18 GPa
- **Fixator**: Variable stiffness (35-140 GPa) for parametric studies

### 2. Finite Element Analysis
- 2D plane stress analysis using 4-node quadrilateral elements
- Adaptive mesh generation with configurable resolution
- Time-stepping simulation (20 steps) modeling healing progression

### 3. Analysis Metrics
- **Von Mises stress** distribution in all material regions
- **Average stress** calculations for bone, callus, and fixator
- **Fracture gap strain** evolution over healing time
- **Stress shielding** effect quantification

### 4. Visualization
- Stress contour plots at each time step
- Comparative stress evolution graphs
- Gap strain progression charts
- Parametric comparison visualizations

## 📈 Output Files

### Generated Results
- `summary_final_results.csv` - Final stress and strain values for all fixator types
- `parametric_comparison.png` - Side-by-side comparison of all simulation results
- Individual simulation folders with:
  - `simulation_results.csv` - Time-series data
  - `stress_shielding.png` - Stress evolution plot
  - `gap_strain.png` - Strain evolution plot
  - `stress_images/` - Stress contour animations

### Key Metrics Tracked
- Final average Von Mises stress in fixator, bone, and callus
- Fracture gap strain evolution
- Stress transfer patterns during healing

## 📊 Results Interpretation

### Stress Shielding Analysis
- **High fixator stiffness** → Increased stress in fixator, reduced bone stress
- **Low fixator stiffness** → More uniform stress distribution, better bone loading

### Gap Strain Evolution
- **Initial high strain** → Gradual reduction as callus stiffens
- **Final strain levels** → Indicator of healing success


## 📚 Documentation

Detailed technical documentation and results are available in the [report](report/) directory, including:
- Complete methodology description
- Numerical results and analysis
- Validation studies
- References to relevant literature
