import numpy as np


Ybus = np.array([[0.04-0.1j, -0.04+0.02j, 0.0+0.08j],
                 [-0.04+0.02j, 0.06-0.1j, 0.0+0.04j],
                 [0.0+0.08j, 0.0+0.04j, 0.03-0.12j]])
Sbus = np.array([100 + 50j, 80 + 30j, 0 + 0j])
V0 = np.array([1.0, 1.0, 1.0])


num_buses = Ybus.shape[0]
V = V0.copy()
P = np.real(Sbus)
Q = np.imag(Sbus)
Pcalc = np.zeros(num_buses)
Qcalc = np.zeros(num_buses)

print(Pcalc)