from neuron import h, gui
import matplotlib.pyplot as plt
from matplotlib import rcParams
from datetime import *
import sys
rcParams['font.family'] = 'SimHei'
plt.rcParams['axes.unicode_minus'] =False
# sys.setrecursionlimit(30000)
current_time = datetime.now()
print(current_time)

h.load_file('stdrun.hoc')
h.cvode_active(1)
cv = h.CVode()
cv.atolscale("DcDt.ng", 1e-22)
cv.atolscale("DcDt.U", 1)
cv.atolscale("DcDt.Z", 1e-6)
cv.jacobian(2)
# cv.maxorder(12)
cv.stiff(0)

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
stim.tdur = 1000
t_vec = h.Vector().record(h._ref_t, 0.1)
c_vec = h.Vector().record(soma(0.5)._ref_cm, 0.1)
z_vec = h.Vector().record(stim._ref_U, 0.1)
v_vec = h.Vector().record(soma(0.5)._ref_v, 0.1)

h.tstop=1500
# h.dt = 1e-6
h.finitialize(-65)
h.run()

end_time = datetime.now()
print(end_time)

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
