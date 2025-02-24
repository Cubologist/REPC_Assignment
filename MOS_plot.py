import numpy as np
import matplotlib.pyplot as plt

# --- Parameters ---
f = 80e3                   # switching frequency [Hz]
T = 1/f                    # switching period [s]
D = 0.7                    # duty cycle
Vin = 45.0                 # input voltage [V]
Vout = 150.0               # output voltage [V]

# Suppose from power-balance you found avg inductor current (boost in CCM):
#   Iout = Vout / R,  Pin = Pout for ideal converter,
#   Iin = (Vout / Vin)*Iout, etc.
# Example: if the load current is 0.375 A, then average inductor current might be ~1.25 A:
I_L_avg = 1.25            # average inductor current [A]

# 20% ripple => peak-to-peak ripple = 0.2 * I_L_avg
delta_IL = 0.2 * I_L_avg   # total p-p ripple
# Let L be whatever value you computed that yields 20% ripple at D=0.7, f=80 kHz, Vin=45V:
#   delta_IL = (Vin * D * T) / L  =>  L = (Vin * D * T) / delta_IL
L = (Vin * D * T) / delta_IL

# Calculate the inductor slopes:
slope_on  =  Vin / L            # when switch is ON (inductor sees +Vin)
slope_off = (Vin - Vout) / L    # when switch is OFF (inductor sees Vin - Vout, negative!)

# Time axis (one switching period)
N  = 500                        # number of points in one period
ts = np.linspace(0, T, N)

# Prepare arrays
vL  = np.zeros_like(ts)
iL  = np.zeros_like(ts)
iSW = np.zeros_like(ts)
iD  = np.zeros_like(ts)

# Build the piecewise waveforms
t_on  = D * T
t_off = T - t_on

# Inductor current starts at (avg - half ripple) at t=0
iL_start = I_L_avg - 0.5*delta_IL
for idx, t in enumerate(ts):
    if t <= t_on:
        # ----- Switch ON -----
        # inductor voltage = Vin
        vL[idx]  = Vin
        # iL(t) rises from iL_start with slope_on
        iL[idx]  = iL_start + slope_on * t
        # switch current = iL(t), diode current = 0
        iSW[idx] = iL[idx]
        iD[idx]  = 0.0
    else:
        # ----- Switch OFF -----
        # inductor voltage = Vin - Vout
        vL[idx]  = Vin - Vout
        # figure out how long we've been in OFF interval
        t_off_elapsed = t - t_on
        # iL at start of OFF = iL(t_on)
        iL_on_to_off = iL_start + slope_on * t_on
        # slope_off applies after t_on
        iL[idx]  = iL_on_to_off + slope_off * t_off_elapsed
        # diode current = inductor current, switch current = 0
        iSW[idx] = 0.0
        iD[idx]  = iL[idx]

# Close the ripple so that iL(end) = iL_start
# (Ideal waveforms should naturally do this if your L matches your ripple spec.)

# --- Plotting ---
fig, axs = plt.subplots(4, 1, figsize=(8,10), sharex=True)

axs[0].plot(ts*1e6, vL, label='v_L(t)')
axs[0].set_ylabel('Inductor Voltage [V]')
axs[0].grid(True)
axs[0].legend()

axs[1].plot(ts*1e6, iL, label='i_L(t)', color='C1')
axs[1].set_ylabel('Inductor Current [A]')
axs[1].grid(True)
axs[1].legend()

axs[2].plot(ts*1e6, iSW, label='i_sw(t)', color='C2')
axs[2].set_ylabel('Switch Current [A]')
axs[2].grid(True)
axs[2].legend()

axs[3].plot(ts*1e6, iD, label='i_d(t)', color='C3')
axs[3].set_ylabel('Diode Current [A]')
axs[3].set_xlabel('Time [µs]')
axs[3].grid(True)
axs[3].legend()

plt.suptitle('Ideal Boost Converter Waveforms in CCM MOS')
plt.tight_layout()
plt.show()
