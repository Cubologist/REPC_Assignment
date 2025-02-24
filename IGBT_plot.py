import numpy as np
import matplotlib.pyplot as plt

# --- Parameters ---
f = 15e3                     # switching frequency [Hz]
T = 1/f                      # switching period [s]
D = 0.85                     # duty cycle
Vout = 300.0                 # output voltage [V]
R = 400.0                    # load resistance [ohms]
Vload_current = Vout / R     # load current [A] = 0.75 A
# Ideal converter power balance => Pin = Pout => Vin * Iin = Vout * Iout
# Also, average inductor current ~ input current in CCM (neglecting losses).

# From Vout = Vin/(1-D), solve for required input voltage:
Vin = Vout * (1 - D)  # ~ 45 V

# Average inductor (input) current = output power / Vin
Iout = Vout / R          # 0.75 A
Pout = Vout * Iout       # 225 W
Iin = Pout / Vin         # ~5 A (average inductor current in CCM)
I_L_avg = Iin            # 5 A

# 25% ripple => peak-to-peak ripple = 0.25 * I_L_avg
delta_IL = 0.25 * I_L_avg   # 1.25 A

# Compute inductance L so that the ripple is as specified:
#   delta_IL = (Vin * D * T) / L
L = (Vin * D * T) / delta_IL

# Slopes for inductor current:
slope_on  =  Vin / L                  # [A/s] when switch is ON
slope_off = (Vin - Vout) / L          # [A/s] when switch is OFF (should be negative)

# Build time axis for one switching period:
N  = 500                             # number of points to plot per period
ts = np.linspace(0, T, N)
t_on  = D * T
t_off = T - t_on

# Arrays for waveforms:
vL  = np.zeros_like(ts)   # Inductor voltage
iL  = np.zeros_like(ts)   # Inductor current
iSW = np.zeros_like(ts)   # Switch (IGBT) current
iD  = np.zeros_like(ts)   # Diode current

# Assume inductor current starts at (I_L_avg - delta_IL/2) at t=0
iL_start = I_L_avg - 0.5 * delta_IL

for idx, t in enumerate(ts):
    if t <= t_on:
        # Switch ON
        vL[idx]  = Vin
        iL[idx]  = iL_start + slope_on * t
        iSW[idx] = iL[idx]
        iD[idx]  = 0.0
    else:
        # Switch OFF
        vL[idx]  = Vin - Vout
        t_off_elapsed = t - t_on
        iL_on_to_off = iL_start + slope_on * t_on
        iL[idx]  = iL_on_to_off + slope_off * t_off_elapsed
        iSW[idx] = 0.0
        iD[idx]  = iL[idx]

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

plt.suptitle('Ideal Boost Converter Waveforms in CCM IGBT')
plt.tight_layout()
plt.show()
