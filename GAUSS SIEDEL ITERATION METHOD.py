import numpy as np

def power_flow_analysis(bus_data, line_data, slack_bus):
    num_buses = len(bus_data)
    num_lines = len(line_data)

    # Create an array to store the bus voltage magnitudes and angles
    voltage = np.ones(num_buses, dtype=complex)
    voltage[slack_bus] = 1.0 + 0.0j  # Slack bus voltage is set to 1.0 per unit magnitude

    # Calculate the admittance matrix
    admittance = np.zeros((num_buses, num_buses), dtype=complex)
    for line in line_data:
        from_bus = line[0]
        to_bus = line[1]
        admittance_val = 1 / (line[2] + line[3] * 1j)
        admittance[from_bus, from_bus] += admittance_val
        admittance[to_bus, to_bus] += admittance_val
        admittance[from_bus, to_bus] -= admittance_val
        admittance[to_bus, from_bus] -= admittance_val

    # Perform Gauss-Seidel iterations
    max_iterations = 100
    epsilon = 1e-6  # Convergence criterion
    for iteration in range(max_iterations):
        max_error = 0.0

        for bus in range(num_buses):
            if bus == slack_bus:
                continue

            # Calculate the bus injection power
            power_injection = bus_data[bus][2] - bus_data[bus][3] * 1j

            # Calculate the bus voltage phasor
            voltage[bus] = power_injection / np.conj(voltage[bus])

            # Calculate the total bus current injection
            total_injection = 0.0j
            for neighbor in range(num_buses):
                total_injection += admittance[bus, neighbor] * voltage[neighbor]

            # Update the bus voltage phasor
            voltage[bus] = power_injection / np.conj(total_injection)

            # Calculate the error in voltage magnitude
            error = abs(voltage[bus] - np.conj(total_injection))
            max_error = max(max_error, error)

        if max_error < epsilon:
            print(f"Converged in {iteration+1} iterations.")
            break

    # Print the bus voltages
    print("Bus Voltages:")
    for bus in range(num_buses):
        magnitude = abs(voltage[bus])
        angle = np.angle(voltage[bus], deg=True)
        print(f"Bus {bus+1}: {magnitude:.4f} ∠ {angle:.2f}°")


# Example data
bus_data = [
    # (bus_number, P_gen, P_load, Q_gen, Q_load)
    (1, 0, 0, 0, 0),
    (2, 2, 0, 1, 0),
    (3, 0, 0.5, 0, 0),
    (4, 1, 0, 0.5, 0),
    (5, 0, 0.2, 0, 0),
]

line_data = [
    # (from_bus, to_bus, resistance, reactance)
    (1, 2, 0.1, 0.2),
    (2, 3, 0.05, 0.15),
    (2, 4, 0.05, 0.1),
    (4, 5, 0.08, 0.25),
]

slack_bus = 0

# Run power flow analysis
power_flow_analysis(bus_data, line_data, slack_bus)

# In this example, bus_data contains the data for each bus in the system, including bus number, active power generation, active power load, reactive power generation, and reactive power load. line_data contains the data for each transmission line, including the "from" bus, "to" bus, resistance, and reactance. The slack_bus variable determines which bus is the slack bus (