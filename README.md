# Gradually Varied Flow (GVF) Solver

A Python-based numerical solver for analyzing gradually varied flow profiles in open channels using the Runge-Kutta 4th order (RK4) method.

## Overview

This tool computes water surface profiles in trapezoidal open channels under gradually varied flow conditions. It automatically determines flow regime (subcritical/supercritical), calculates critical and normal depths, and simulates the water surface profile using numerical integration.

## Features

- **Automatic Flow Classification**: Identifies mild (M), steep (S), or critical (C) slope profiles
- **Numerical Depth Calculation**: Solves for normal depth (yn) and critical depth (yc) using Newton-Raphson iteration
- **RK4 Integration**: High-accuracy 4th-order Runge-Kutta method for profile computation
- **Intelligent Direction Logic**: Automatically determines upstream/downstream integration based on flow regime
- **Visual Output**: Generates publication-quality plots with reference depths and transition zones
- **Robust Error Handling**: Prevents singularities near critical depth and dry bed conditions

## Theory

The solver integrates the gradually varied flow equation:

```
dy/dx = (S₀ - Sf) / (1 - Fr²)
```

Where:
- `S₀` = bed slope
- `Sf` = friction slope (from Manning's equation)
- `Fr` = Froude number

### Flow Regimes

| Regime | Condition | Integration Direction | Control Location |
|--------|-----------|----------------------|------------------|
| Subcritical | y > yc | Upstream (negative dx) | Downstream |
| Supercritical | y < yc | Downstream (positive dx) | Upstream |

## Installation

### Requirements

```bash
pip install numpy matplotlib
```

### Supported Python Versions
- Python 3.6+

## Usage

### Basic Example

Run the script and follow the interactive prompts:

```bash
python gvf_solver.py
```

### Input Parameters

The program will request the following inputs:

1. **Channel Parameters**:
   - `Q`: Discharge (m³/s)
   - `b`: Bottom width (m)
   - `m`: Side slope (horizontal:vertical ratio)
   - `n`: Manning's roughness coefficient
   - `S₀`: Bed slope

2. **Boundary Conditions**:
   - `y_start`: Starting depth (m)
   - `x_start`: Starting position (m)
   - `sim_length`: Length to simulate (m)

### Example Session

```
--- GRADUALLY VARIED FLOW SOLVER (RK4) ---

[ CHANNEL PARAMETERS ]
Discharge Q (m^3/s) [e.g., 20]: 20
Bottom Width b (m) [e.g., 10]: 10
Side Slope m (H:V) [e.g., 1.5]: 1.5
Manning's n [e.g., 0.015]: 0.015
Bed Slope S0 [e.g., 0.0005]: 0.0005

[ CALCULATED REFERENCE DEPTHS ]
Normal Depth (yn)   : 1.2345 m
Critical Depth (yc) : 0.9876 m
Slope Type: MILD (M-profile)

[ BOUNDARY CONDITION ]
Enter Starting Depth y (m): 1.5
Enter Starting Position x (m): 0
Enter Length to Simulate (m): 1000

-> Subcritical Flow detected (y > yc).
-> Integrating BACKWARDS (Upstream).
```

## Output

The program generates a matplotlib plot showing:

- **Blue solid line**: Computed water surface profile
- **Green dashed line**: Normal depth (yn)
- **Red dash-dot line**: Critical depth (yc)
- **Yellow shaded area**: Transition zone between yn and yc
- **Black line**: Channel bottom

## Profile Classifications

### Mild Slope (yn > yc)
- **M1**: Depth decreases downstream (backwater curve)
- **M2**: Depth increases downstream (drawdown curve)
- **M3**: Supercritical flow above normal depth

### Steep Slope (yn < yc)
- **S1**: Subcritical flow on steep slope
- **S2**: Rapid transition to normal depth
- **S3**: Supercritical acceleration

## Technical Details

### Numerical Methods

1. **Newton-Raphson Iteration**: Used to solve implicit equations for normal and critical depths with tolerance of 1×10⁻⁶

2. **RK4 Integration**: Fourth-order accurate spatial integration with adaptive stopping conditions

3. **Singularity Handling**: Automatically detects and stops computation near critical depth (|1 - Fr²| < 0.01)

### Channel Geometry

The solver supports trapezoidal channels with the following formulations:

- **Area**: `A = (b + m·y)·y`
- **Wetted Perimeter**: `P = b + 2y√(1 + m²)`
- **Top Width**: `T = b + 2m·y`
- **Hydraulic Radius**: `R = A/P`

## Limitations

- Assumes steady, gradually varied flow
- Does not model hydraulic jumps or rapids
- Requires uniform channel geometry
- Step size fixed at 5m (can be modified in code)
- Trapezoidal cross-sections only

## Troubleshooting

**Issue**: "Normal depth did not converge"
- **Solution**: Adjust initial guess in `solve_normal_depth()` or check input parameters

**Issue**: "Simulation stopped early"
- **Solution**: Flow approached critical depth or boundary. This is expected behavior near transitions.

**Issue**: Unrealistic profiles
- **Solution**: Verify Manning's n values and ensure boundary condition (y_start) is physically reasonable

---

**Note**: Always verify computational results against field measurements or established hydraulic design software for critical applications.
