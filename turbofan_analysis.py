"""
Ideal Turbofan Engine Thermodynamic Analysis
Calculates cycle states, thrust, and efficiencies, and generates a T-s diagram.
"""

import math
import numpy as np
import matplotlib.pyplot as plt

def get_constants():
    """Returns thermodynamic constants and ambient conditions at 41,000 ft."""
    return {
        "gamma_air": 1.4,
        "R_air": 287.0,
        "Cp_air": 1005.0,
        "gamma_gas": 1.33,
        "R_gas": 287.0,
        "Cp_gas": 1150.0,
        "p_0": 17870.0,             # Ambient pressure (Pa)
        "T_0": -56.5 + 273.15,      # Ambient temperature (K)
        "Ma_cr": 0.94               # Cruise Mach number
    }

def get_engine_parameters():
    """Returns the main design parameters of the turbofan engine."""
    return {
        "bypass_ratio": 5.0,
        "overall_pr": 23.0,
        "mass_flow_rate": 110.0,    # Total mass flow rate (kg/s)
        "fan_pr": 1.6,              # Fan pressure ratio
        "T_turbine_inlet": 1994.0   # Turbine inlet temperature (K)
    }

def calculate_cycle(c, params):
    """Performs the thermodynamic calculations for the ideal turbofan cycle."""
    results = {}
    
    # Ambient conditions and velocities
    a_0 = math.sqrt(c["gamma_air"] * c["R_air"] * c["T_0"])
    V_cr = c["Ma_cr"] * a_0
    results["V_cr"] = V_cr

    # Free stream stagnation (Total) properties
    p_t_0 = c["p_0"] * pow((1 + (c["gamma_air"] - 1) / 2 * c["Ma_cr"]**2), (c["gamma_air"] / (c["gamma_air"] - 1)))
    T_t_0 = c["T_0"] * (1 + (c["gamma_air"] - 1) / 2 * c["Ma_cr"]**2)

    # Mass flows
    m_core = params["mass_flow_rate"] / (1 + params["bypass_ratio"])
    m_bypass = params["mass_flow_rate"] - m_core

    # Fan & Bypass Stream (Isentropic)
    p_t_1 = p_t_0 * params["fan_pr"]
    T_t_1 = T_t_0 * pow(params["fan_pr"], (c["gamma_air"] - 1) / c["gamma_air"])
    
    p_bypass_exit = c["p_0"]
    Ma_bypass_exit = math.sqrt((2 / (c["gamma_air"] - 1)) * (pow(p_t_1 / p_bypass_exit, (c["gamma_air"] - 1) / c["gamma_air"]) - 1))
    T_bypass_exit = T_t_1 * pow(p_bypass_exit / p_t_1, (c["gamma_air"] - 1) / c["gamma_air"])
    V_bypass_exit = Ma_bypass_exit * math.sqrt(c["gamma_air"] * c["R_air"] * T_bypass_exit)

    # Core Compressor (Isentropic)
    compressor_pr = params["overall_pr"] / params["fan_pr"]
    p_t_2 = p_t_1 * compressor_pr
    T_t_2 = T_t_1 * pow(compressor_pr, (c["gamma_air"] - 1) / c["gamma_air"])
    
    # Combustor (Isobaric)
    p_t_3 = p_t_2
    T_t_3 = params["T_turbine_inlet"]
    Q_combustion = m_core * (c["Cp_gas"] * T_t_3 - c["Cp_air"] * T_t_2)

    # Turbine (Work balance: W_turbine = W_compressor + W_fan)
    W_compressor = m_core * c["Cp_air"] * (T_t_2 - T_t_1)
    W_fan = params["mass_flow_rate"] * c["Cp_air"] * (T_t_1 - T_t_0)
    W_turbine = W_compressor + W_fan
    
    T_t_4 = T_t_3 - W_turbine / (m_core * c["Cp_gas"])
    p_t_4 = p_t_3 * pow(T_t_4 / T_t_3, c["gamma_gas"] / (c["gamma_gas"] - 1))

    # Core Nozzle (Isentropic)
    p_core_exit = c["p_0"]
    Ma_core_exit = math.sqrt((2 / (c["gamma_gas"] - 1)) * (pow(p_t_4 / p_core_exit, (c["gamma_gas"] - 1) / c["gamma_gas"]) - 1))
    T_core_exit = T_t_4 * pow(p_core_exit / p_t_4, (c["gamma_gas"] - 1) / c["gamma_gas"])
    V_core_exit = Ma_core_exit * math.sqrt(c["gamma_gas"] * c["R_gas"] * T_core_exit)

    # Thrust & Efficiencies
    thrust_bypass = m_bypass * (V_bypass_exit - V_cr)
    thrust_core = m_core * (V_core_exit - V_cr)
    total_thrust = thrust_bypass + thrust_core

    delta_KE = 0.5 * m_bypass * (V_bypass_exit**2 - V_cr**2) + 0.5 * m_core * (V_core_exit**2 - V_cr**2)
    propulsive_power = total_thrust * V_cr
    
    thermal_eff = delta_KE / Q_combustion
    propulsive_eff = propulsive_power / delta_KE
    overall_eff = thermal_eff * propulsive_eff

    # Package results for printing and plotting
    results.update({
        "m_core": m_core, "m_bypass": m_bypass, "p_t_0": p_t_0, "T_t_0": T_t_0,
        "p_t_1": p_t_1, "T_t_1": T_t_1, "p_t_2": p_t_2, "T_t_2": T_t_2,
        "p_t_3": p_t_3, "T_t_3": T_t_3, "p_t_4": p_t_4, "T_t_4": T_t_4,
        "Q_combustion": Q_combustion, "V_bypass_exit": V_bypass_exit, "V_core_exit": V_core_exit,
        "total_thrust": total_thrust, "thrust_bypass": thrust_bypass, "thrust_core": thrust_core,
        "thermal_eff": thermal_eff, "propulsive_eff": propulsive_eff, "overall_eff": overall_eff,
        "T_bypass_exit": T_bypass_exit, "T_core_exit": T_core_exit, "T_0": c["T_0"]
    })
    
    return results

def print_performance_report(res, params):
    """Prints the calculated engine performance metrics."""
    print("=" * 50)
    print(" IDEAL TURBOFAN ENGINE PERFORMANCE REPORT")
    print("=" * 50)
    print(f"Cruise Velocity        : {res['V_cr']:.2f} m/s")
    print(f"Core Mass Flow Rate    : {res['m_core']:.2f} kg/s")
    print(f"Bypass Mass Flow Rate  : {res['m_bypass']:.2f} kg/s")
    print("-" * 50)
    print(f"Total Thrust           : {res['total_thrust']:.2f} N")
    print(f"Bypass Thrust          : {res['thrust_bypass']:.2f} N")
    print(f"Core Thrust            : {res['thrust_core']:.2f} N")
    print(f"Specific Thrust        : {res['total_thrust'] / params['mass_flow_rate']:.2f} N/(kg/s)")
    print("-" * 50)
    print(f"Thermal Efficiency     : {res['thermal_eff'] * 100:.2f} %")
    print(f"Propulsive Efficiency  : {res['propulsive_eff'] * 100:.2f} %")
    print(f"Overall Efficiency     : {res['overall_eff'] * 100:.2f} %")
    print("=" * 50)

def plot_ts_diagram(res, c):
    """Generates and displays the Temperature-Entropy (T-s) diagram."""
    # Assuming ambient entropy s_0 = 0 for reference
    s_0 = 0
    s_t_0 = s_0    # Isentropic ram
    s_t_1 = s_t_0  # Isentropic fan
    s_t_2 = s_t_1  # Isentropic compressor

    # Combustor curve (Isobaric heat addition)
    T_combustor = np.linspace(res['T_t_2'], res['T_t_3'], 100)
    s_combustor = s_t_2 + c["Cp_gas"] * np.log(T_combustor / res['T_t_2'])
    s_t_3 = s_combustor[-1]

    s_t_4 = s_t_3  # Isentropic turbine
    s_core_exit = s_t_4 
    s_bypass_exit = s_t_1 

    # Plot coordinates
    T_core_pts = [res['T_0'], res['T_t_0'], res['T_t_1'], res['T_t_2'], res['T_t_3'], res['T_t_4'], res['T_core_exit']]
    s_core_pts = [s_0, s_t_0, s_t_1, s_t_2, s_t_3, s_t_4, s_core_exit]

    T_bypass_pts = [res['T_t_1'], res['T_bypass_exit']]
    s_bypass_pts = [s_t_1, s_bypass_exit]

    plt.figure(figsize=(10, 7))
    
    # Plot cycles
    plt.plot(s_core_pts[:4], T_core_pts[:4], 'b-', linewidth=2.5, label='Compression (Isentropic)')
    plt.plot(s_combustor, T_combustor, 'r-', linewidth=2.5, label='Combustion (Isobaric)')
    plt.plot(s_core_pts[4:], T_core_pts[4:], 'g-', linewidth=2.5, label='Core Expansion (Isentropic)')
    plt.plot(s_bypass_pts, T_bypass_pts, 'c--', linewidth=2, label='Bypass Expansion')

    # Scatter station points
    plt.scatter(s_core_pts, T_core_pts, color='black', zorder=5)

    # Annotations
    offsets = {
        s_0: (res['T_0'] - 50, '0 (Ambient)'),
        s_t_0: (res['T_t_0'], 't0 (Inlet)'),
        s_t_1: (res['T_t_1'] + 30, 't1 (Fan Exit)'),
        s_t_2: (res['T_t_2'], 't2 (Compressor Exit)'),
        s_t_3: (res['T_t_3'] + 30, 't3 (Turbine Inlet)'),
        s_t_4: (res['T_t_4'], 't4 (Turbine Exit)'),
        s_core_exit: (res['T_core_exit'] - 50, 'Core Exit')
    }
    
    for s_val, (T_val, label) in offsets.items():
        plt.text(s_val, T_val, f' {label}', fontsize=10, fontweight='bold')

    plt.title('Ideal Turbofan Engine: T-s Diagram', fontsize=14, fontweight='bold')
    plt.xlabel('Entropy [s] (J/kg*K)', fontsize=12)
    plt.ylabel('Total Temperature [T] (K)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper left', fontsize=11)
    
    plt.tight_layout()
    plt.show()

def main():
    constants = get_constants()
    parameters = get_engine_parameters()
    
    results = calculate_cycle(constants, parameters)
    print_performance_report(results, parameters)
    plot_ts_diagram(results, constants)

if __name__ == "__main__":
    main()
