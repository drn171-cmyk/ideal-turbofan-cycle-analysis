# Ideal Turbofan Engine Thermodynamic Analysis

This repository contains a Python-based numerical solver for the thermodynamic cycle analysis of an ideal turbofan engine. It evaluates the performance of the engine at cruise altitude and visualizes the complete Temperature-Entropy (T-s) diagram.

## Key Features
* **Thermodynamic Cycle Solver:** Calculates stagnation temperatures and pressures across all engine stations (Fan, Compressor, Combustor, Turbine, Nozzles).
* **Performance Metrics:** Computes Core & Bypass Thrust, Specific Fuel Consumption (SFC), Thermal Efficiency, Propulsive Efficiency, and Overall Efficiency.
* **Visualization:** Generates an accurate T-s diagram highlighting isentropic compression/expansion and isobaric heat addition.

## Engineering Parameters Investigated
* **Flight Altitude:** 41,000 ft
* **Cruise Mach Number:** 0.94
* **Bypass Ratio (BPR):** 5.0
* **Overall Pressure Ratio (OPR):** 23.0
* **Turbine Inlet Temperature (TIT):** 1994 K

## How to Run
1. Install dependencies: `pip install numpy matplotlib`
2. Execute the script: `python turbofan_analysis.py`

## Developer
* **Diren Gürgül** - Mechanical Engineering Student at Istanbul Technical University (ITU)
