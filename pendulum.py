# Packages 
import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt 

# System Parameters - from user guide
g = 9.81        # m/s**2 - gravity 
Rm = 8.4        # ohms - terminal resistance
Kt = 0.042      # N-m/A - torque const
Km = 0.042      # V/rad/s - motor back emf const
mr = 0.095      # kg - rotary arm mass
Lr = 0.085      # m - rotary arm length
mp = 0.024      # kg - pendulum mass
Lp = 0.129      # m - pendulum length

# Moment of inertia - rod about end
Jp = (1/3)*mp*Lp*Lp
Jr = (1/3)*mr*Lr*Lr

# Damping - idk put a function here later, maybe move into dynamics if it depends on current val
Dr = 0
Dp = 0

# Time step 
h = 0.025       # seconds

# System dynamics 
def pendulum_dynamics(current_state, input_voltage):
    """
    System equations for dynamics of rotary pendulum
    """
    x1, x2, x3, x4 = current_state
    Vm = input_voltage

    x1_dot = x3
    x2_dot = x4
    x3_dot = -(2*(8*Dr*Jp*Rm*x3 - 2*Kt*Lp**2*Vm*mp - 8*Jp*Kt*Vm + 8*Jp*Km*Kt*x3 + 2*Dr*Lp**2*Rm*mp*x3 + 2*Km*Kt*Lp**2*mp*x3 + Lp**3*Lr*Rm*mp**2*x4**2*sin(x2) + 4*Dp*Lp*Lr*Rm*mp*x4*cos(x2) - Lp**3*Lr*Rm*mp**2*x3**2*cos(x2)**2*sin(x2) - 2*Lp**2*Lr*Rm*g*mp**2*cos(x2)*sin(x2) + 4*Jp*Lp*Lr*Rm*mp*x4**2*sin(x2) + Lp**4*Rm*mp**2*x3*x4*cos(x2)*sin(x2) + 4*Jp*Lp**2*Rm*mp*x3*x4*cos(x2)*sin(x2)))/(Rm*(Lp**4*mp**2 + 16*Jp*Jr - Lp**4*mp**2*cos(x2)**2 + 4*Lp**2*Lr**2*mp**2 + 4*Jp*Lp**2*mp + 16*Jp*Lr**2*mp + 4*Jr*Lp**2*mp - 4*Lp**2*Lr**2*mp**2*cos(x2)**2 - 4*Jp*Lp**2*mp*cos(x2)**2))
    x4_dot = -(16*Dp*Jr*Rm*x4 - 2*Lp**3*Rm*g*mp**2*sin(x2) + 4*Dp*Lp**2*Rm*mp*x4 + 16*Dp*Lr**2*Rm*mp*x4 - 8*Jr*Lp*Rm*g*mp*sin(x2) + Lp**4*Rm*mp**2*x3**2*cos(x2)**3*sin(x2) + 2*Lp**3*Rm*g*mp**2*cos(x2)**2*sin(x2) - 4*Dp*Lp**2*Rm*mp*x4*cos(x2)**2 - Lp**4*Rm*mp**2*x3**2*cos(x2)*sin(x2) - 8*Lp*Lr**2*Rm*g*mp**2*sin(x2) - 8*Kt*Lp*Lr*Vm*mp*cos(x2) + 8*Dr*Lp*Lr*Rm*mp*x3*cos(x2) + 8*Km*Kt*Lp*Lr*mp*x3*cos(x2) - 4*Lp**2*Lr**2*Rm*mp**2*x3**2*cos(x2)*sin(x2) + 4*Lp**2*Lr**2*Rm*mp**2*x4**2*cos(x2)*sin(x2) - 4*Jr*Lp**2*Rm*mp*x3**2*cos(x2)*sin(x2) + 4*Lp**3*Lr*Rm*mp**2*x3*x4*cos(x2)**2*sin(x2))/(Rm*(Lp**4*mp**2 + 16*Jp*Jr - Lp**4*mp**2*cos(x2)**2 + 4*Lp**2*Lr**2*mp**2 + 4*Jp*Lp**2*mp + 16*Jp*Lr**2*mp + 4*Jr*Lp**2*mp - 4*Lp**2*Lr**2*mp**2*cos(x2)**2 - 4*Jp*Lp**2*mp*cos(x2)**2))
 
    new_state = np.array([x1_dot, x2_dot, x3_dot, x4_dot])

    return new_state

# Runge-kutta integrate 

def rk4_integrate(x, u, h):
    """
    Runge–Kutta method of integration
    """

    k1 = pendulum_dynamics(x, u)
    k2 = pendulum_dynamics(x + h*k1/2, u)
    k3 = pendulum_dynamics(x + h*k2/2, u)
    k4 = pendulum_dynamics(x + h*k3, u)

    rk4_sol = x + (h/6)*(k1+2*k2+2*k3+k4)

    return rk4_sol

# Basic simulation i guess hahahheh

T = 50 #seconds - length of simulation 

initial_state = np.array([0, 0.1, 0, 0]) #starting at just theta != 0 gives nothing lmao 

def simple_simulation(initial_state, input_u, time_length, time_step):
    """
    Really basic simulation
    """
    state = initial_state
    #something to store values in 
    sim_data = np.append(np.insert(state, 0, 0), input_u)

    for i in range(int(time_length/h)):
        state = rk4_integrate(state, input_u, time_step)   #run 1 interation 

        # store sim data - probably a better way but this is easy
        state_val = np.append(np.insert(state, 0, i*time_step), input_u)

        sim_data = np.append(sim_data, state_val)
    

    #rehape array to be pretty 
    sim_data = np.reshape(sim_data, (-1, 6))

    # make a plot 
    plt.plot(sim_data[:, 0], sim_data[:, 1])
    plt.plot(sim_data[:, 0], sim_data[:, 2])
    plt.plot(sim_data[:, 0], sim_data[:, 3])
    plt.plot(sim_data[:, 0], sim_data[:, 4])
    plt.legend((r'$\theta$', r'$\alpha$', r'$\dot{\theta}$', r'$\dot{\alpha}$'),
           loc='upper center', shadow=True)
    plt.show()

    plt.plot(sim_data[:, 0], sim_data[:, 1])
    plt.hlines(np.pi, 0, T, linestyles='dashed')
    plt.hlines(-np.pi, 0, T, linestyles='dashed')
    plt.show()



simple_simulation(initial_state, 0, T, h)

