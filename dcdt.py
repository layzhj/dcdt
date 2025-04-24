from neuron import h
import matplotlib.pyplot as plt
from matplotlib import rcParams
import time

rcParams['font.family'] = 'SimHei'
plt.rcParams['axes.unicode_minus'] =False

start_time = time.time()


h.load_file('stdgui.hoc')
h.cvode_active(1)
cv = h.CVode()
cv.atolscale("DcDt.ng", 1e-22)
cv.atolscale("DcDt.U", 1)
cv.atolscale("DcDt.Z", 1e-6)


A = 300e3
f = 500
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
stim.rel_Zmin = -0.1

stim.tbegin = 100
stim.tdur = 100
t_vec = h.Vector().record(h._ref_t)
c_vec = h.Vector().record(soma(0.5)._ref_cm)
z_vec = h.Vector().record(stim._ref_Z)
ng_vec = h.Vector().record(stim._ref_ng)
q_vec = h.Vector().record(stim._ref_q)
stm = h.Vector().record(stim._ref_stm)

h.tstop=300
# h.dt = 1e-6
h.finitialize(-65)
try:
    h.run()
except:
    pass

end_time = time.time()
print('运行时长：', end_time - start_time)

t_array = t_vec.as_numpy()
z_array = z_vec.as_numpy()
c_array = c_vec.as_numpy()
ng_array = ng_vec.as_numpy()
q_array = q_vec.as_numpy()
stm_array = stm.as_numpy()
print(z_array)


fig, ax = plt.subplots(5, 1)

ax[0].plot(t_array, z_array, color='orange')
ax[0].set_ylabel('Z (um)')

ax[1].plot(t_array, ng_array, color='orange')
ax[1].set_ylabel('ng (mol)')

ax[2].plot(t_array, q_array, color='orange')
ax[2].set_ylabel('Charge ($nC/cm^{2}$)')

ax[3].plot(t_array, stm_array, color='orange')
ax[3].set_ylabel('Driver (Pa)')

ax[4].plot(t_array, c_array, color='orange', label='Capacitance')
ax[4].set_xlabel('Time (ms)')
ax[4].set_ylabel('Capacitance ($uF/cm^{2}$)')
plt.show()
