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
theta0 = 45*(np.pi/180) # initial angle (radians)
omega0 = 0.0        # initial angular velocity (rad/s)  
y0 = [theta0, omega0]
t_span = (0, 10)  # time span for the simulation (seconds)
t_eval = np.linspace(t_span[0], t_span[1], 1000)
# Solve the ODE
sol_linear = solve_ivp(pendulum_ode_linear, t_span, y0, args=(g, l), t_eval=t_eval)


# Plot the results
plt.figure(figsize=(10, 5)) 
plt.plot(sol_linear.t, sol_linear.y[0]*180/np.pi, label='Linear')
plt.xlabel('Time (s)')
plt.ylabel('Angle (deg)')


## add nonlinear solution
# def pendulum_ode_nonlinear(t, y, g, l):
#     theta, omega = y
#     dtheta_dt = omega
#     domega_dt = -(g / l) * np.sin(theta)
#     return [dtheta_dt, domega_dt]
# sol_nonlinear = solve_ivp(pendulum_ode_nonlinear, t_span, y0, args=(g, l), t_eval=t_eval)
# ## add nonlinear solution
# plt.plot(
#     sol_nonlinear.t, sol_nonlinear.y[0] * 180 / np.pi, label="Nonlinear", linestyle="--"
# )
# plt.legend()


# using a different solver
# sol_radau = solve_ivp(pendulum_ode_nonlinear, t_span, y0, args=(g, l), t_eval=t_eval, method="Radau")
# plt.plot(
#     sol_radau.t, sol_radau.y[0] * 180 / np.pi, label="Nonlinear Radau", linestyle=":"
# )
# plt.legend()


# zoom in x-axis
# plt.xlim(t_span[1]-10, t_span[1])