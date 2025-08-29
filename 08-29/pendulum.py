## numerically solve a simple linearized pendulum
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def pendulum_ode_linear(t, y, g, l):
    theta, omega = y
    dtheta_dt = omega
    domega_dt = -(g / l) * theta
    return [dtheta_dt, domega_dt]

# Parameters
g = 9.81  # acceleration due to gravity (m/s^2)
l = 1.0   # length of the pendulum (m)
theta0 = np.pi / 4  # initial angle (radians)
omega0 = 0.0        # initial angular velocity (rad/s)  
y0 = [theta0, omega0]
t_span = (0, 10)  # time span for the simulation (seconds)
t_eval = np.linspace(t_span[0], t_span[1], 100)
# Solve the ODE
sol = solve_ivp(pendulum_ode_linear, t_span, y0, args=(g, l), t_eval=t_eval)

# Plot the results
plt.figure(figsize=(10, 5)) 
plt.plot(sol.t, sol.y[0], label='Angle (rad)')

