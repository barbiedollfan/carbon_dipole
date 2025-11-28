import math
import pandas as pd
import matplotlib.pyplot as plt

a = 0.0245  # in pm^-1 
r_e = 116.0 # put your value here (in pm) 
D_e = 1326.8 # D_e = . . . you'll need to figure that out. It needs to be in zJ = 10⁻²¹ J 
m_O = 26.57     # mass of oxygen atom in yg 
m_C = 19.93     # mass of carbon atom in yg 

# Equilibrium positions 

x1_eq = -r_e # oxygen 1 sitting at x = -r_e 
x2_eq = 0 # carbon sitting at x = 0 
x3_eq = +r_e # oxygen 2 sitting at x = +r_e 

TIME_STEP = 0.0001     # this worked well for me 

def update_x(x, v, a): 
    new_x = x + v*TIME_STEP + 0.5*a*TIME_STEP**2 
    return new_x 

def update_v(v, a): 
    new_v = v + a*TIME_STEP 
    return new_v 

def compute_dipole(x_1, x_2, x_3): 
    d = -0.65*x_1 + 1.30*x_2 -0.65*x_3 
    return d 

def F(r):
    exp = math.exp(-a * (r - r_e))
    return -2.0 * a * D_e * (1.0 - exp) * exp

records_1 = []  # records motion of x_1 
records_2 = []  # records motion of x_2 
records_3 = []  # records motion of x_3 
records_d = []  # records dipole moment 


# Initial positions (pm) 

delta = 1.51 #(estimated from energy absorbed from a single photon) 

x_1 = x1_eq + delta 
x_2 = x2_eq - (2*m_O/m_C) * delta 
x_3 = x3_eq + delta 

# Initial velocities (pm/fs = 10³ m/s) 
v_1 = 0.0 
v_2 = 0.0 
v_3 = 0.0 

# Initial dipole moment 
d =  compute_dipole(x_1, x_2, x_3) 
t = 0.0 
count = 0.0

T_MAX = 200     # a good starting value, but experiment 

while t < T_MAX: 
    # compute F12 
    r12 = abs(x_2 - x_1) 
    F12 = F(r12) 

    # compute F23 
    r23 = abs(x_3 - x_2) 
    F23 = F(r23) 

    # compute accel of mass 1 (oxygen 1) 
    a_1 = -F12 / m_O    # Fnet = ma for oxygen #1 

    # compute accel of mass 2 (carbon) 
    a_2 = (F12 - F23) / m_C # Fnet = ma for carbon 

    # compute a3 
    a_3 = F23 / m_O     # Fnet = ma for oxygen #2 

    # Store position, velocity, acceleration data of first Oxygen 
    records_1.append( 
        {'Time': t, 'Position':x_1, 'Velocity': v_1, 'Acceleration': a_1} 
    ) 

    # Store position, velocity, acceleration data of Carbon 
    records_2.append( 
        {'Time': t, 'Position':x_2, 'Velocity': v_2, 'Acceleration': a_2} 
    ) 

    # Store position, velocity, acceleration data of second Oxygen 
    records_3.append( 
        {'Time': t, 'Position':x_3, 'Velocity': v_3, 'Acceleration': a_3} 
    ) 

    # Store dipole moment 
    records_d.append( 
        {'Time': t, 'DipoleMoment': d} 
    ) 

    # Update x_1 
    x_1 = update_x(x_1, v_1, a_1) 

    # Update x_2 
    x_2 = update_x(x_2, v_2, a_2) 

    # Update x_3 
    x_3 = update_x (x_3, v_3, a_3) 

    # Update velocity 1 
    v_1 = update_v(v_1, a_1) 

    # Update velocity 2 
    v_2 = update_v(v_2, a_2) 

    # Update velocity 3 
    v_3 = update_v(v_3, a_3) 

    # compute new dipole moment 
    d =  compute_dipole(x_1, x_2, x_3) 

    # Update time 
    t += TIME_STEP 
    count += 1 

print(f"Done after {count} iterations.") 

# Unpack the records into pandas dataframes 
x_1_data = pd.DataFrame.from_records(records_1) 
x_2_data = pd.DataFrame.from_records(records_2) 
x_3_data = pd.DataFrame.from_records(records_3) 

dipole_data = pd.DataFrame.from_records(records_d) 

# Common axis label 
TIME_LABEL = r"Time (fs = $10^{-15}$ s)" 

def plot_time_series(df, y, title, ylabel, fig_num): 
    """Generic helper to plot a time series from a DataFrame.""" 
    plt.figure(fig_num) 
    df.plot(x="Time", y=y) 
    plt.title(title) 
    plt.xlabel(TIME_LABEL) 
    plt.ylabel(ylabel) 
    plt.tight_layout() 
    plt.savefig(f"{fig_num}.png")
    plt.show() 

# Species and their dataframes 
species = [ 
    ("Oxygen 1", x_1_data), 
    ("Carbon",   x_2_data), 
    ("Oxygen 2", x_3_data), 
] 

# Quantities to plot (y-column name, y-axis label) 
quantities = [ 
    ("Position",     r"Position (pm = $10^{-12}$ m)"), 
    ("Velocity",     r"Velocity (pm/fs = $10^3$ m/s)"), 
    ("Acceleration", r"Acceleration ($\mathrm{pm}/\mathrm{fs}^2 = 10^{18}$ $\mathrm{m}/\mathrm{s}^2$)"), 
] 

# Generate figures 1–9 for position/velocity/acceleration 
fig_num = 1 

for y, ylabel in quantities: 
    for title, df in species: 
        plot_time_series(df, y=y, title=title, ylabel=ylabel, fig_num=fig_num) 
        fig_num += 1 

# Plot dipole moment as a function of time (figure 10) 
plt.figure(fig_num) 
dipole_data.plot(x="Time", y="DipoleMoment") 
plt.title("Dipole Moment") 
plt.xlabel(TIME_LABEL) 
plt.ylabel(r"Dipole Moment (e·pm = $1.602\times 10^{-31}$ C·m)") 
plt.tight_layout() 
plt.savefig("dipole.png")
plt.show() 
