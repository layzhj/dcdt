from neuron import h
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'SimHei'
plt.rcParams['axes.unicode_minus'] =False


h.load_file('stdgui.hoc')
h.cvode_active(1)
cv = h.CVode()

# 

A = 300e3
f = 500e3
soma = h.Section(name='soma')
soma.L = 10
soma.diam = 10
soma.cm = 1
ref_cm = soma(0.5)._ref_cm
stim = h.DcDt(soma(0.5))
stim._ref_c = ref_cm
stim.A = A
stim.f = f
stim.Delta = 1.2735619230798087e-03
stim.z0 = 6.3168827639957465e-6


stim.tbegin = 100
stim.tdur = 100
t_vec = h.Vector().record(h._ref_t)
c_vec = h.Vector().record(soma(0.5)._ref_cm)
z_vec = h.Vector().record(stim._ref_U)
v_vec = h.Vector().record(soma(0.5)._ref_v)

h.tstop=300
# h.dt = 0.025
h.finitialize(-65)
h.run()

t_array = t_vec.as_numpy()
z_array = z_vec.as_numpy()
c_array = c_vec.as_numpy()
v_array = v_vec.as_numpy()
print(z_array)

fig, ax = plt.subplots(3, 1)

ax[0].plot(t_array, z_array, color='orange')
ax[0].set_ylabel('Z (um)')

ax[1].plot(t_array, v_array, color='orange')
ax[1].set_ylabel('V (mV)')

ax[2].plot(t_array, c_array, color='orange', label='Capacitance')
ax[2].set_xlabel('Time (ms)')
ax[2].set_ylabel('Capacitance ($uF/cm^{2}$)')
plt.show()