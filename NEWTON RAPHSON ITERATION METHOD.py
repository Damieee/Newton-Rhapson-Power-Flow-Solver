# The Newton-Raphson method is an iterative technique widely used for load flow analysis in power systems. It converges to a solution by iteratively linearizing and solving a set of nonlinear power flow equations. Here's a simplified explanation of the Newton-Raphson load flow iteration method:

# Start with an initial guess for the voltage magnitudes and angles at all buses.
# Compute the complex power injections at each bus based on the current voltage estimates and system parameters.
# Linearize the power flow equations by computing the Jacobian matrix.
# Solve the linearized equations to obtain the correction values for voltage angles and magnitudes.
# Update the voltage estimates by adding the correction values.
# Repeat steps 2 to 5 until convergence is achieved.
# Convergence is typically determined by comparing the mismatch between the calculated and specified power injections against a tolerance level.
# Here's a high-level Python implementation of the Newton-Raphson load flow iteration method:

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

def newton_raphson_load_flow(Ybus, Sbus, V0, max_iter=20, tol=1e-6):
    # Initialize variables
    num_buses = Ybus.shape[0]
    V = V0.copy().astype(complex)  # Convert V0 to complex data type
    P = np.real(Sbus)
    Q = np.imag(Sbus)

    # Perform Newton-Raphson load flow iterations
    for iteration in range(max_iter):
        # Compute complex power injections
        S = V * np.conj(np.dot(Ybus, V))

        # Compute Jacobian matrix
        dS_dVm = np.diag(V) * np.conj(Ybus * np.conj(V))
        dS_dVa = 1j * np.diag(V) * np.conj(Ybus * V)
        J11 = np.real(dS_dVa) - np.imag(dS_dVm)
        J12 = np.real(dS_dVm) + np.imag(dS_dVa)
        J21 = np.imag(dS_dVa) + np.real(dS_dVm)
        J22 = np.imag(dS_dVm) - np.real(dS_dVa)
        J = np.vstack((np.hstack((J11, J12)), np.hstack((J21, J22))))

        # Compute power mismatch
        mismatch = np.concatenate((P - np.real(S), Q - np.imag(S)))

        # Solve linearized equations
        dx = spsolve(csc_matrix(J), mismatch)

        # Update voltage angles and magnitudes
        dtheta = dx[:num_buses]
        dV = dx[num_buses:]
        V += dV
        V *= np.exp(1j * dtheta)

        # Check convergence
        if np.max(np.abs(mismatch)) < tol:
            break

    VoltageLevel_VoltageAngle=[]
    for data in V:
        # Separating real and imaginary parts
        real_part = np.real(data)
        # Calculate the angle in degrees
        angle_deg = np.angle(data, deg=True)
        VoltageLevel_VoltageAngle.append(f"{[real_part, angle_deg]}")


    return VoltageLevel_VoltageAngle



# Example usage
Ybus = np.array([[0.04-0.1j, -0.04+0.02j, 0.0+0.08j],
                 [-0.04+0.02j, 0.06-0.1j, 0.0+0.04j],
                 [0.0+0.08j, 0.0+0.04j, 0.03-0.12j]])
Sbus = np.array([100 + 50j, 80 + 30j, 0 + 0j])
V0 = np.array([1.0, 1.0, 1.0])
V = newton_raphson_load_flow(Ybus, Sbus, V0)
print("Voltage Magnitudes and Angles:")
print(V)

# In this implementation, the newton_raphson_load_flow function takes the admittance matrix Ybus, complex power injections Sbus, and initial voltage estimates V0 as inputs. It iteratively performs the Newton-Raphson load flow method until convergence is achieved or the maximum number of iterations is reached. The resulting voltage magnitudes and angles are then returned.

# The example usage demonstrates how to call the newton_raphson_load_flow function with a sample admittance matrix Ybus, complex power injections Sbus, and initial voltage estimates V0. The resulting voltage magnitudes and angles are then printed.

# Please note that this is a simplified implementation for educational purposes, and it may not handle certain cases or system configurations. It's always recommended to use validated and tested power flow libraries or consult domain experts for critical power flow analysis tasks.