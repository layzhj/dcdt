from neuron import h
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'SimHei'
plt.rcParams['axes.unicode_minus'] =False


h.load_file('stdgui.hoc')
h.cvode_active(1)
cv = h.CVode()

A = 500e6
f = 300e3
soma = h.Section(name='soma')
soma.L = 10
soma.diam = 10
soma.cm = 1
ref_cm = soma(0.5)._ref_cm
stim = h.DcDt(soma(0.5))
stim._ref_c = ref_cm
stim.A = A
stim.f = f

stim.tbegin = 100
stim.tdur = 100
t_vec = h.Vector().record(h._ref_t)
c_vec = h.Vector().record(soma(0.5)._ref_cm)
z_vec = h.Vector().record(stim._ref_Z)

h.tstop=300
# h.dt = 0.025
h.finitialize(-65)
h.run()

t_array = t_vec.as_numpy()
z_array = z_vec.as_numpy()
c_array = c_vec.as_numpy()
print(z_array)


plt.plot(t_array, c_array, color='orange', label='Capacitance')
plt.xlabel('Time (ms)')
plt.ylabel('Capacitance ($uF/cm^{2}$)')
plt.show()