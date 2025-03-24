# Power-System-State-Estimation-with-WLS
This repository contains a MATLAB implementation of a power system state estimation algorithm using the Weighted Least Squares (WLS) method. The code is designed to work with IEEE standard benchmark systems (e.g., IEEE 14-bus, 30-bus) via MATPOWER case files, enabling analysis of voltage magnitudes and angles under noisy measurement conditions. It includes tools for measurement generation, noise simulation, state estimation, and bad data detection using a chi-squared (χ²) test.

## Features
- **State Estimation**: Estimates bus voltage magnitudes and angles using WLS based on noisy measurements.
- **Measurement Generation**: Extracts voltage magnitudes, power injections, and power flows from MATPOWER case files with customizable bus and branch selections.
- **Noise Simulation**: Adds realistic multiplicative Gaussian noise to true measurements, reflecting typical power system uncertainties.
- **Bad Data Detection**: Implements a χ²-test to identify anomalies in measurement data.
- **Visualization**: Plots estimated vs. true voltage magnitudes and angles, and χ²-test results over time.

## Requirements
- **MATLAB**: Tested with MATLAB R2023a (earlier versions may work but are untested).
- **MATPOWER**: Required for loading IEEE benchmark case files (e.g., `case14.m`). Install via [MATPOWER’s website](https://matpower.org/) or [GitHub](https://github.com/MATPOWER/matpower).
- **Dependencies**:
  - `loadcase` (from MATPOWER)
  - Standard MATLAB functions (e.g., `chi2inv`, `randn`, `plot`)
 
## Installation
1. **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/Power-System-State-Estimation-with-WLS.git
   cd Power-System-State-Estimation-with-WLS
2. **Install MATPOWER**:
   Download MATPOWER from [here](https://matpower.org/download/) and add it to your MATLAB path:
    
   ```matlab
   addpath('/path/to/matpower');

3. **Verify Setup**:
   Run `test_matpower` in MATLAB to ensure MATPOWER is correctly installed.

## Files
- `computeBBus.m`: Computes the B-bus (Shunt Admittance Matrix).
- `computeYBus.m`: compute the bus admittance matrix (Ybus).
- `extractLineData.m`: Extracts line datas.
- `extractZData.m`: Generates measurement structure.
- `powerflow.m`: true measurements for the WLS estimation.
- `wls.m`: Implements the WLS state estimation algorithm.
- `chi_squared_test.m`: Performs χ²-test for bad data detection.
- `plot_results.m`: Plots estimated vs. true voltage magnitudes and angles.

## Methodology
- **Measurement Uncertainty**: Uses realistic noise levels (0.2% for voltages, 1% for injections, 0.5% for flows) based on modern devices (e.g., PMUs, SCADA).
- **WLS Algorithm**: Iteratively solves for system state using a Jacobian matrix and measurement residuals.
- **Bad Data Detection**: Compares the WLS objective function (`J`) against a χ² threshold.

## Contributing
Contributions are welcome! Please:
1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/your-feature`).
3. Commit changes (`git commit -m 'Add your feature'`).
4. Push to the branch (`git push origin feature/your-feature`).
5. Open a pull request.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
