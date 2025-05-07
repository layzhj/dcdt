from neuron import h
import matplotlib.pyplot as plt
from matplotlib import rcParams
import time
import pandas as pd

rcParams['font.family'] = 'SimHei'
plt.rcParams['axes.unicode_minus'] =False

start_time = time.time()

save_path = './hoc_long.csv'

h.load_file('stdgui.hoc')
h.cvode_active(1)
cv = h.CVode()
cv.atolscale("DcDt.ng", 1e-22)
cv.atolscale("DcDt.U", 1)
cv.atolscale("DcDt.Z", 1e-6)
# cv.atol(1e-10)


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
stim.rel_Zmin = -0.49

stim.tbegin = 100
stim.tdur = 50
# t_vec = h.Vector().record(h._ref_t)
# c_vec = h.Vector().record(soma(0.5)._ref_cm)
# z_vec = h.Vector().record(stim._ref_Z)
# ng_vec = h.Vector().record(stim._ref_ng)
# q_vec = h.Vector().record(stim._ref_q)
# stm = h.Vector().record(stim._ref_stm)
# U_vec = h.Vector().record(stim._ref_U)
t_vec = h.Vector().record(h._ref_t, 0.1)
c_vec = h.Vector().record(soma(0.5)._ref_cm, 0.1)
z_vec = h.Vector().record(stim._ref_Z, 0.1)
ng_vec = h.Vector().record(stim._ref_ng, 0.1)
q_vec = h.Vector().record(stim._ref_q, 0.1)
stm = h.Vector().record(stim._ref_stm, 0.1)
U_vec = h.Vector().record(stim._ref_U, 0.1)

h.tstop=150
# h.dt = 1e-6
h.finitialize(-65)
try:
    h.run()
except:
    pass

end_time = time.time()
print('continue:', end_time - start_time)

t_array = t_vec.as_numpy()
z_array = z_vec.as_numpy()
c_array = c_vec.as_numpy()
ng_array = ng_vec.as_numpy()
q_array = q_vec.as_numpy()
stm_array = stm.as_numpy()
U_array = U_vec.as_numpy()
print(z_array)

# df = pd.DataFrame({'t': t_array, 'z': z_array, 'ng': ng_array, 'q': q_array, 'stm': stm_array, 'U': U_array, 'c': c_array})
# df.to_csv(save_path, index=False)

fig, ax = plt.subplots(6, 1)

ax[0].plot(t_array, z_array, color='orange')
ax[0].set_ylabel('Z (um)')

ax[1].plot(t_array, ng_array, color='orange')
ax[1].set_ylabel('ng (mol)')

ax[2].plot(t_array, q_array, color='orange')
ax[2].set_ylabel('Charge ($nC/cm^{2}$)')

ax[3].plot(t_array, stm_array, color='orange')
ax[3].set_ylabel('Driver (Pa)')

ax[4].plot(t_array, U_array, color='orange')
ax[4].set_ylabel('dZ/dt (um/ms)')

ax[5].plot(t_array, c_array, color='orange', label='Capacitance')
ax[5].set_xlabel('Time (ms)')
ax[5].set_ylabel('Capacitance ($uF/cm^{2}$)')
plt.show()
