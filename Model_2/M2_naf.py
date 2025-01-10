# Script for model 2 of the falling cat - Peano's model (Naf)

import numpy as np
import matplotlib.pyplot as plt

class PeanoModel:
    def __init__(self, segments=2, inertia=None, mass=5.0, gravity=9.81, lever_arm=0.1, initial_height=10.0):
        self.segments = segments  # Number of body segments (e.g., front and back)
        self.segment_angles = np.zeros(segments, dtype=float)  # Angles of each segment
        self.segment_angular_velocities = np.zeros(segments, dtype=float)  # Angular velocities of each segment
        self.inertia = np.array(inertia if inertia is not None else np.ones(segments), dtype=float)  # Moments of inertia
        self.mass = mass  # Mass of the cat (kg)
        self.gravity = gravity  # Acceleration due to gravity (m/s^2)
        self.lever_arm = lever_arm  # Distance to the center of mass (m)
        self.height = initial_height  # Initial height (m)
        self.vertical_velocity = 0.0  # Initial vertical velocity (m/s)

    def compute_torque(self, angle):
        """
        Compute the torque acting on a segment due to gravity.
        :param angle: Angle of the segment (radians).
        :return: Torque (Nm).
        """
        return -self.mass * self.gravity * self.lever_arm * np.sin(angle)

    def update_angles(self, delta_time):
        """
        Update segment angles considering torque and angular momentum conservation.
        :param delta_time: Time increment for the simulation.
        """
        # Compute torques for each segment
        torques = np.array([self.compute_torque(angle) for angle in self.segment_angles])

        # Compute angular accelerations
        angular_accelerations = torques / self.inertia

        # Update angular velocities
        self.segment_angular_velocities += angular_accelerations * delta_time

        # Conserve angular momentum
        total_angular_momentum = np.sum(self.inertia * self.segment_angular_velocities)
        correction_factor = total_angular_momentum / np.sum(self.inertia)
        self.segment_angular_velocities -= correction_factor / self.inertia

        # Update angles based on corrected angular velocities
        self.segment_angles += self.segment_angular_velocities * delta_time

    def update_height(self, delta_time):
        """
        Update the vertical position of the cat.
        :param delta_time: Time increment for the simulation.
        """
        self.vertical_velocity += self.gravity * delta_time
        self.height -= self.vertical_velocity * delta_time

    def get_segment_angles(self):
        """Return the current angles of all segments."""
        return self.segment_angles

    def has_hit_ground(self):
        """Check if the cat has hit the ground."""
        return self.height <= 0

# Simulation parameters
segments = 2
inertia = [1, 1]  # Equal moments of inertia for simplicity
mass = 5.0  # Average mass of a cat (kg)
gravitational_acceleration = 9.81  # Gravity (m/s^2)
lever_arm = 0.1  # Distance to the center of mass (m)
initial_height = 100 # Initial height (m)
initial_angular_velocities = [1, -1]  # Initial velocities ensuring zero total angular momentum
simulation_time = 100 # Maximum simulation time in seconds
delta_time = 0.01  # Time step

def run_simulation():
    peano_model = PeanoModel(segments=segments, inertia=inertia, mass=mass, gravity=gravitational_acceleration, lever_arm=lever_arm, initial_height=initial_height)
    peano_model.segment_angular_velocities = np.array(initial_angular_velocities, dtype=float)

    times = [0]
    segment_angles_over_time = []
    heights_over_time = [initial_height]

    print("Starting simulation...")

    while times[-1] < simulation_time:
        # Update dynamics
        peano_model.update_angles(delta_time)
        peano_model.update_height(delta_time)

        # Record data
        times.append(times[-1] + delta_time)
        segment_angles_over_time.append(peano_model.get_segment_angles().copy())
        heights_over_time.append(peano_model.height)

        # Check if the cat has hit the ground
        if peano_model.has_hit_ground():
            print(f"Cat hit the ground at time: {times[-1]:.2f} seconds")
            break

    # Convert data to arrays for plotting
    times = np.array(times)
    segment_angles_over_time = np.array(segment_angles_over_time)
    heights_over_time = np.array(heights_over_time)

    print("Simulation completed. Plotting results...")

    # Plot results
    plt.figure(figsize=(8, 6))
    for segment in range(segments):
        plt.plot(times[:-1], segment_angles_over_time[:, segment], label=f"Segment {segment + 1} Angle")
    plt.axhline(y=0, color='k', linestyle='--', label="Zero Rotation")
    plt.xlabel("Time (s)")
    plt.ylabel("Segment Angles (Radians)")
    plt.title("Peano Model Simulation: Angles During Fall")
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.plot(times, heights_over_time, label="Height (m)")
    plt.axhline(y=0, color='r', linestyle='--', label="Ground Level")
    plt.xlabel("Time (s)")
    plt.ylabel("Height (m)")
    plt.title("Peano Model Simulation: Height During Fall")
    plt.legend()
    plt.grid()
    plt.show()

run_simulation()