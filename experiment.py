#packages 
import numpy as np
from numpy import sin, cos 

state_values = {
    'x1': 2,    # theta
    'x2': 3,    # alpha
    'x3': 3,    # theta velocity
    'x4': 5,    # alpha velocity
}

global_parameters = {
    "mr": 0.095,            # kg - rotary arm mass
    "Lr": 0.085,            # m - rotary arm length
    "mp": 0.024,            # kg - pendulum mass
    "Lp": 0.129,            # m - pendulum length
    "Jr": 11111111,
    "Jp": 1111111,
    "Kt": 0.042,            # N-m/A - torque const
    "Km": 0.042,            # V/rad/s - motor back emf const
    "Rm": 8.4,              # ohms - terminal resistance 
    "g": 9.81,              # m/s**2 - gravity 
}


# System Parameters - from user guide
g = 9.81        # m/s**2 - gravity 
Rm = 8.4        # ohms - terminal resistance
Kt = 0.042      # N-m/A - torque const
Km = 0.042      # V/rad/s - motor back emf const
mr = 0.095      # kg - rotary arm mass
Lr = 0.085      # m - rotary arm length
mp = 0.024      # kg - pendulum mass
Lp = 0.129      # m - pendulum length

#moment of inertia - rod about end
Jp = (1/3)*mp*Lp*Lp
Jr = (1/3)*mr*Lr*Lr

#Damping
Dr = 0
Dp = 0

def pendulum_dynamics(current_state, input_voltage):
    """
    System equations for dynamics of pendulum
    """

    x1, x2, x3, x4 = current_state
    Vm = input_voltage

    x1_dot = x3
    x2_dot = x4
    x3_dot = -(2*(8*Dr*Jp*Rm*x3 - 2*Kt*Lp**2*Vm*mp - 8*Jp*Kt*Vm + 8*Jp*Km*Kt*x3 + 2*Dr*Lp**2*Rm*mp*x3 + 2*Km*Kt*Lp**2*mp*x3 + Lp**3*Lr*Rm*mp**2*x4**2*sin(x2) + 4*Dp*Lp*Lr*Rm*mp*x4*cos(x2) - Lp**3*Lr*Rm*mp**2*x3**2*cos(x2)**2*sin(x2) - 2*Lp**2*Lr*Rm*g*mp**2*cos(x2)*sin(x2) + 4*Jp*Lp*Lr*Rm*mp*x4**2*sin(x2) + Lp**4*Rm*mp**2*x3*x4*cos(x2)*sin(x2) + 4*Jp*Lp**2*Rm*mp*x3*x4*cos(x2)*sin(x2)))/(Rm*(Lp**4*mp**2 + 16*Jp*Jr - Lp**4*mp**2*cos(x2)**2 + 4*Lp**2*Lr**2*mp**2 + 4*Jp*Lp**2*mp + 16*Jp*Lr**2*mp + 4*Jr*Lp**2*mp - 4*Lp**2*Lr**2*mp**2*cos(x2)**2 - 4*Jp*Lp**2*mp*cos(x2)**2))
    x4_dot = -(16*Dp*Jr*Rm*x4 - 2*Lp**3*Rm*g*mp**2*sin(x2) + 4*Dp*Lp**2*Rm*mp*x4 + 16*Dp*Lr**2*Rm*mp*x4 - 8*Jr*Lp*Rm*g*mp*sin(x2) + Lp**4*Rm*mp**2*x3**2*cos(x2)**3*sin(x2) + 2*Lp**3*Rm*g*mp**2*cos(x2)**2*sin(x2) - 4*Dp*Lp**2*Rm*mp*x4*cos(x2)**2 - Lp**4*Rm*mp**2*x3**2*cos(x2)*sin(x2) - 8*Lp*Lr**2*Rm*g*mp**2*sin(x2) - 8*Kt*Lp*Lr*Vm*mp*cos(x2) + 8*Dr*Lp*Lr*Rm*mp*x3*cos(x2) + 8*Km*Kt*Lp*Lr*mp*x3*cos(x2) - 4*Lp**2*Lr**2*Rm*mp**2*x3**2*cos(x2)*sin(x2) + 4*Lp**2*Lr**2*Rm*mp**2*x4**2*cos(x2)*sin(x2) - 4*Jr*Lp**2*Rm*mp*x3**2*cos(x2)*sin(x2) + 4*Lp**3*Lr*Rm*mp**2*x3*x4*cos(x2)**2*sin(x2))/(Rm*(Lp**4*mp**2 + 16*Jp*Jr - Lp**4*mp**2*cos(x2)**2 + 4*Lp**2*Lr**2*mp**2 + 4*Jp*Lp**2*mp + 16*Jp*Lr**2*mp + 4*Jr*Lp**2*mp - 4*Lp**2*Lr**2*mp**2*cos(x2)**2 - 4*Jp*Lp**2*mp*cos(x2)**2))
 
    new_state = [x1_dot, x2_dot, x3_dot, x4_dot]

    return new_state


pendulum_dynamics([1, 3, 4, 5], 2)


def linearised_system_dyamics(parameters, current_state):

    x1, x2, x3, x4 = current_state

    J_T = parameters["Jp"]*parameters["mp"]*parameters["Lr"]**2 + parameters["Jr"]*parameters["Jp"] + 0.25*parameters["Jr"]*parameters["mp"]*parameters["Lp"]**2

    print("J", J_T)

    x1_dot = x3
    x2_dot = x4
    x3_dot = 6
    x4_dot = 7

    theta_ddot = (1/J_T)*(parameters["Jp"]+0.25*parameters["mp"]*parameters["Lp"]**2)*(torque-parameters["Dr"]*x3)+0.5*(1/J_T)*parameters["mp"]*parameters["Lp"]*parameters["Lr"]*(0.5*parameters["mp"]*parameters["Lp"]*parameters["g"])*x2-parameters["Dp"]*x4

    return J_T

linearised_system_dyamics(global_parameters)