# Fracture Fixation FEA Simulation

A comprehensive finite element analysis (FEA) toolkit for investigating the biomechanical behavior of bone-fixator systems during fracture healing. This project models stress transfer mechanisms and analyzes how different fixator stiffnesses influence the healing process through computational simulation.

## 🎯 Project Overview

This project leverages 2D finite element analysis to explore the critical biomechanical phenomenon of "stress transfer" during the bone fracture healing process. The simulation models:

-   **A Bone-Fixator System:** Incorporating time-dependent callus maturation to simulate the healing progression.
-   **Stress Redistribution:** Analyzing how mechanical stresses are distributed among the bone, fixator, and developing callus tissues.
-   **Fracture Gap Strain Evolution:** Tracking changes in strain within the fracture gap under various fixator rigidities.
-   **Parametric Studies:** Conducting comparative analyses using flexible, standard, and rigid fixator configurations.

## 🏗️ Project Structure

The project is organized as follows:

```txt
.
├── main.py                      # Main script to run simulations and generate results
├── README.md                    # This README file
├── utils/                       # Directory for core simulation modules
│   ├── fea_core.py              # Core finite element analysis logic
│   ├── materials.py             # Definitions of material properties
│   ├── simulation.py            # Simulation workflow and control
│   ├── analysis.py              # Post-processing: stress and strain calculations
│   └── plot_utils.py            # Utilities for generating plots and visualizations
├── output_advanced/             # Default directory for generated simulation results
│   ├── Flexible/                # Results for the flexible fixator
│   ├── Standard/                # Results for the standard fixator
│   └── Rigid/                   # Results for the rigid fixator
├── assets/                      # Contains model structure figure and generated animations
└── report/                      # LaTeX source documentation
```

## 🚀 Quick Start

### Prerequisites

Ensure you have Python 3 installed. Then, install the necessary libraries:

```bash
pip install numpy matplotlib scipy pandas
```

### Running the Simulation

Execute the main script from the project's root directory:

```bash
python main.py
```

## 📊 Key Features

### 1. Material Modeling

The simulation incorporates distinct material properties:
-   **Bone:** Modeled as an elastic material with an Elastic Modulus of 18 GPa and a Poisson's ratio of 0.3.
-   **Callus:** Simulates the healing process with a time-dependent stiffening behavior, where its Elastic Modulus gradually increases from an initial 0.1 GPa up to 18 GPa, matching mature bone.
-   **Fixator:** Subject to parametric study with three configurations based on Elastic Modulus:
    -   **Flexible Fixator:** 35 GPa (35.0e9 Pa)
    -   **Standard Fixator:** 70 GPa (70.0e9 Pa)
    -   **Rigid Fixator:** 140 GPa (140.0e9 Pa)

### 2. Finite Element Analysis

-   **Methodology:** Employs 2D plane stress analysis.
-   **Elements:** Utilizes 4-node quadrilateral elements for meshing the domain.
-   **Meshing:** Features adaptive mesh generation with configurable resolution to balance accuracy and computational cost.
-   **Simulation Steps:** A time-stepping simulation (default: 20 steps) models the progression of fracture healing and callus maturation.

### 3. Analysis Metrics

The simulation tracks and calculates several key biomechanical indicators:
-   **Von Mises Stress:** Distribution across all material regions (bone, callus, fixator).
-   **Average Stress:** Mean stress values calculated for bone, callus, and fixator components.
-   **Fracture Gap Strain:** Evolution of strain within the fracture gap over the healing period.
-   **Stress Shielding Effect:** Quantified by comparing stress levels in the bone under different fixator stiffnesses.

### 4. Visualization

Comprehensive visualization tools are integrated to aid in the interpretation of results:
-   **Stress Contour Plots:** Generated for each time step, illustrating stress distribution.
-   **Comparative Stress Evolution Graphs:** Plotting average stress in different components over time for each fixator type.
-   **Gap Strain Progression Charts:** Showing how fracture gap strain changes throughout the healing simulation.
-   **Parametric Comparison Visualizations:** Side-by-side plots for easy comparison of outcomes from different fixator stiffnesses.

## 📈 Output Files

Upon successful execution, `main.py` generates the following in the `output_advanced/` directory:

-   **`summary_final_results.csv`**: A CSV file summarizing the final average Von Mises stress in the fixator, bone, and callus, along with the final fracture gap strain for each fixator type.
-   **`parametric_comparison.png`**: An image file providing a side-by-side visual comparison of key results across all simulated fixator types.
-   **Individual Simulation Folders (`Flexible/`, `Standard/`, `Rigid/`)**:
    -   `simulation_results.csv`: Time-series data for average stresses and gap strain for that specific simulation run.
    -   `stress_shielding.png`: A plot showing the evolution of average stress in bone, callus, and fixator over time.
    -   `gap_strain.png`: A plot illustrating the evolution of fracture gap strain over time.
    -   `stress_images/`: A directory containing individual frames of stress contour plots, which can be compiled into an animation.

## 🔬 Results Interpretation

The simulation outputs provide insights into:

### Stress Shielding Analysis

-   **High Fixator Stiffness (Rigid Fixator):** Tends to bear a larger portion of the mechanical load, leading to increased stress within the fixator itself and potentially reducing the mechanical stimulus (stress) experienced by the bone. This is known as stress shielding.
-   **Low Fixator Stiffness (Flexible Fixator):** Allows for more load sharing with the bone and callus, resulting in a more uniform stress distribution. This can promote better bone loading, which is often considered beneficial for healing.

### Gap Strain Evolution

-   **Initial Phase:** The fracture gap typically experiences high strain levels at the beginning of the healing process when the callus is soft.
-   **Healing Progression:** As the callus matures and stiffens, the strain in the fracture gap gradually reduces.
-   **Final Strain Levels:** The magnitude of strain at the end of the simulation can be an indicator of the mechanical stability achieved and potentially correlates with healing success. Interfragmentary strain theory suggests optimal ranges of strain for different stages of healing.

## 📚 Documentation

Detailed technical documentation, theoretical background, and in-depth analysis of results are available in the [report/](report/) directory, typically as a LaTeX-generated PDF document.