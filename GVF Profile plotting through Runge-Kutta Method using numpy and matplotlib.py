import numpy as np
import matplotlib.pyplot as plt
import sys

# ==========================================
# 1. GEOMETRY & HYDRAULIC CLASSES
# ==========================================

class Channel:
    def __init__(self, Q, n, S0, b, m):
        self.Q = Q      # Discharge (m^3/s)
        self.n = n      # Manning's roughness
        self.S0 = S0    # Bed slope
        self.b = b      # Bottom width (m)
        self.m = m      # Side slope (H:V)
        self.g = 9.81   # Gravity

    def area(self, y):
        return (self.b + self.m * y) * y

    def perimeter(self, y):
        return self.b + 2 * y * np.sqrt(1 + self.m**2)

    def top_width(self, y):
        return self.b + 2 * self.m * y

    def hydraulic_radius(self, y):
        A = self.area(y)
        P = self.perimeter(y)
        return A / P if P > 0 else 0

# ==========================================
# 2. SOLVING NORMAL & CRITICAL DEPTH
# ==========================================

def solve_normal_depth(ch):
    """
    Solves Manning's Equation for yn using Newton-Raphson.
    Target: Q - (1/n)*A*R^(2/3)*S0^(1/2) = 0
    """
    y = 4.0  # Initial guess
    tol = 1e-6
    max_iter = 100

    for _ in range(max_iter):
        A = ch.area(y)
        P = ch.perimeter(y)
        R = A / P
        T = ch.top_width(y)
        
        # Function value (Manning Discrepancy)
        Q_calc = (1.0 / ch.n) * A * (R**(2/3)) * (ch.S0**0.5)
        f = Q_calc - ch.Q
        
        # Derivative (finite difference approximation for simplicity)
        dy_small = 1e-5
        A_d = ch.area(y+dy_small)
        P_d = ch.perimeter(y+dy_small)
        R_d = A_d / P_d
        Q_calc_d = (1.0 / ch.n) * A_d * (R_d**(2/3)) * (ch.S0**0.5)
        df = (Q_calc_d - Q_calc) / dy_small
        
        if abs(f) < tol:
            return y
        
        y = y - f / df
        if y <= 0: y = 0.1 # Prevent negative depth during iteration
        
    print("Warning: Normal depth did not converge.")
    return y

def solve_critical_depth(ch):
    """
    Solves Froude=1 for yc using Newton-Raphson.
    Target: (Q^2 * T) / (g * A^3) - 1 = 0
    """
    y = 2.0 # Initial guess
    tol = 1e-6
    max_iter = 100
    
    for _ in range(max_iter):
        A = ch.area(y)
        T = ch.top_width(y)
        
        if A <= 0: A = 0.1
        
        # Froude Number Squared
        Fr2 = (ch.Q**2 * T) / (ch.g * A**3)
        f = Fr2 - 1.0
        
        # Derivative (finite difference)
        dy_small = 1e-5
        A_d = ch.area(y+dy_small)
        T_d = ch.top_width(y+dy_small)
        Fr2_d = (ch.Q**2 * T_d) / (ch.g * A_d**3)
        df = (Fr2_d - Fr2) / dy_small
        
        if abs(f) < tol:
            return y
            
        y = y - f / df
        if y <= 0: y = 0.1

    print("Warning: Critical depth did not converge.")
    return y

# ==========================================
# 3. RUNGE-KUTTA (RK4) IMPLEMENTATION
# ==========================================

def get_dy_dx(x, y, ch):
    """
    The Differential Equation: dy/dx = (S0 - Sf) / (1 - Fr^2)
    """
    if y <= 0.05: return 0 # Min depth cap

    A = ch.area(y)
    P = ch.perimeter(y)
    R = A / P
    T = ch.top_width(y)

    # 1. Friction Slope (Sf)
    # Sf = (n^2 * Q^2) / (A^2 * R^(4/3))
    Sf = (ch.n**2 * ch.Q**2) / (A**2 * np.power(R, 4.0/3.0))

    # 2. Froude Number (Fr)
    Fr2 = (ch.Q**2 * T) / (ch.g * A**3)

    # 3. Singularity Check
    if abs(1 - Fr2) < 0.01:
        # We are hitting critical depth (denominator -> 0).
        # We stop the slope calculation to prevent infinity.
        return 0 

    dydx = (ch.S0 - Sf) / (1 - Fr2)
    return dydx

def solve_profile(ch, x_start, y_start, length, step_size):
    x_vals = [x_start]
    y_vals = [y_start]
    
    current_x = x_start
    current_y = y_start
    
    # Determine number of steps
    # We use abs() because step_size might be negative
    steps = int(length / abs(step_size))
    
    print(f"Starting Simulation: x={x_start}, y={y_start}, dx={step_size}")
    
    for _ in range(steps):
        # RK4 Constants
        k1 = get_dy_dx(current_x, current_y, ch)
        k2 = get_dy_dx(current_x + step_size/2, current_y + k1*step_size/2, ch)
        k3 = get_dy_dx(current_x + step_size/2, current_y + k2*step_size/2, ch)
        k4 = get_dy_dx(current_x + step_size, current_y + k3*step_size, ch)
        
        next_y = current_y + (step_size / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        
        # Stop if depth becomes unrealistic or hits critical depth
        if next_y <= 0 or abs(get_dy_dx(current_x, next_y, ch)) == 0:
            print(f"Simulation stopped early at x={current_x:.2f} (Approached Critical Depth or Dry Bed)")
            break

        current_x += step_size
        current_y = next_y
        
        x_vals.append(current_x)
        y_vals.append(current_y)
        
    return x_vals, y_vals

# ==========================================
# 4. MAIN USER INTERFACE
# ==========================================

def main():
    print("--- GRADUALLY VARIED FLOW SOLVER (RK4) ---")
    
    # --- A. User Inputs ---
    try:
        print("\n[ CHANNEL PARAMETERS ]")
        Q = float(input("Discharge Q (m^3/s) [e.g., 20]: ") or 20)
        b = float(input("Bottom Width b (m) [e.g., 10]: ") or 10)
        m = float(input("Side Slope m (H:V) [e.g., 1.5]: ") or 1.5)
        n = float(input("Manning's n [e.g., 0.015]: ") or 0.015)
        S0 = float(input("Bed Slope S0 [e.g., 0.0005]: ") or 0.0005)
        
        ch = Channel(Q, n, S0, b, m)
        
        # --- B. Calculate Reference Depths ---
        yn = solve_normal_depth(ch)
        yc = solve_critical_depth(ch)
        
        print(f"\n[ CALCULATED REFERENCE DEPTHS ]")
        print(f"Normal Depth (yn)   : {yn:.4f} m")
        print(f"Critical Depth (yc) : {yc:.4f} m")
        
        if yn > yc:
            print("Slope Type: MILD (M-profile)")
        elif yn < yc:
            print("Slope Type: STEEP (S-profile)")
        else:
            print("Slope Type: CRITICAL (C-profile)")
            
        # --- C. Boundary Conditions ---
        print("\n[ BOUNDARY CONDITION ]")
        y_start = float(input(f"Enter Starting Depth y (m): "))
        x_start = float(input(f"Enter Starting Position x (m): "))
        sim_length = float(input(f"Enter Length to Simulate (m): "))
        
        # --- D. Automatic Direction Logic ---
        # If depth > critical, it is subcritical -> control is downstream -> calc upstream (Negative dx)
        # If depth < critical, it is supercritical -> control is upstream -> calc downstream (Positive dx)
        
        dx = 5.0 # default magnitude
        
        if y_start > yc:
            print("\n-> Subcritical Flow detected (y > yc).")
            print("-> Integrating BACKWARDS (Upstream).")
            dx = -abs(dx) # Ensure negative
        else:
            print("\n-> Supercritical Flow detected (y < yc).")
            print("-> Integrating FORWARDS (Downstream).")
            dx = abs(dx) # Ensure positive

        # --- E. Run Simulation ---
        x_plot, y_plot = solve_profile(ch, x_start, y_start, sim_length, dx)
        
        # --- F. Plotting ---
        plt.figure(figsize=(10, 6))
        
        # 1. Plot Water Surface
        plt.plot(x_plot, y_plot, 'b-', linewidth=2, label='Water Surface Profile')
        
        # 2. Plot Reference Lines
        plt.axhline(yn, color='green', linestyle='--', label=f'Normal Depth ({yn:.2f}m)')
        plt.axhline(yc, color='red', linestyle='-.', label=f'Critical Depth ({yc:.2f}m)')
        
        # 3. Plot Bed Slope (Relative)
        # We plot the bed relative to the water depth. 
        # For visualization, we treat the bed as z=0 in the plot, 
        # but to show the slope effect we can calculate bed elevation Z.
        
        # Create Bed Elevation Array
        # Assuming the datum is at x=0 or x=L depending on flow
        # Simple approach: Plot Water Level (Z + y) vs Distance
        
        # Let's stick to Specific Energy style: Depth vs Distance 
        # (Standard for profile classification)
        plt.axhline(0, color='k', linewidth=3, label='Channel Bottom')

        # Fill zones
        plt.fill_between([min(x_plot), max(x_plot)], yn, yc, color='yellow', alpha=0.1, label='Transition Zone')
        
        plt.title(f"GVF Profile | Q={Q}, n={n}, S0={S0}")
        plt.xlabel("Distance along Channel (m)")
        plt.ylabel("Flow Depth (m)")
        plt.legend()
        plt.grid(True, alpha=0.5)
        plt.tight_layout()
        plt.show()
        
    except ValueError:
        print("Error: Please enter valid numeric values.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()