# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 07:17:26 2021

@author: joshu
"""

from hhg2 import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.collections import LineCollection
from scipy.signal import find_peaks
import cmath
from tkinter import *
import json
from tkinter.filedialog import asksaveasfile
from scipy.signal import savgol_filter
#%%
I_p = 21.5645 * e / E_h
omega_0s = 2e6 * np.pi * hbar * c / (E_h * lambda_0)
avo = 6.02214e23

dE = 2e6 * np.pi * c / lambda_0 * (hbar / e)
E_0 = 0.16

Omega_cut = np.arange(30, 801, 20)
sigma_cut = [2.6462e5,2.2118e5,1.8494e5,1.4062e5,1.0346e5,7.4288e4,5.4935e4,4.2866e4,\
             3.2727e4,2.6314e4,2.12e4,1.7354e4,1.4426e4,1.2079e4,1.0202e4,8707,\
             7502,6517,5695,4986,4394,3894,3467,3102,2788,2516,2279,2071,1890,\
             1730,1588,1458,1339,1241,1156,1076,1003,936.7,875.1]
smigma_cut = np.asarray(sigma_cut) * 1e-4 * 20.18 / (a_0**2 * avo)
sigma_cut = np.asarray(sigma_cut) * 1e-4 * 20.18 / avo
lambda_0 = np.arange(0.1, 2.6, 0.1)
lambda_0 = lambda_0 * 1e-6

plt.plot(Omega_cut,smigma_cut)
plt.yscale('log')
plt.show()


pts = []

pts = []
z = []

for i in lambda_0:
    etas = []
    for j,k in enumerate(sigma_cut):
        omega_0 = 2 * np.pi * c / i 
        x = HHG(10, 0.16*E_h/(e*a_0), omega_0, 100, 17, 1e30, k, 21.5645, 5*i/c)
        eta = x.eta_cutoff()
        etas.append(eta)
        pt = [i,Omega_cut[j],eta]
        pts.append(pt)
    z.append(etas)

#x = hhg.HHG(10, E_0, omega_0, L, harmonic, n_density, sigma, I_p, T)
#%%
from hhg import *
import numpy as np
import matplotlib.pyplot as plt

I_p = 21.5645 * e / E_h
avo = 6.02214e23
E_0 = 0.16

Omega_cut = np.arange(30, 801, 20)
sigma_cut = [8.3523e5,3.1238e5,1.4701e5,8.3218e4,5.0042e4,3.1533e4,2.0657e4, \
             1.4235e4,1.0386e4,7753,5944,4658,3613,2897,2359,1945,1623,1367,1162, \
             992.2,853.4,739.1,644.3,565,498.1,441.4,393.1,350.5,314.2,282.8, \
             255.4,230.3,207.7,189.5,173.8,159.2,146.2,134.4,123.8]
smigma_cut = np.asarray(sigma_cut) * 1e-4 * 20.18 / (a_0**2 * avo)
sigma_cut = np.asarray(sigma_cut) * 1e-4 * 20.18 / avo
lambda_0 = np.arange(0.1, 2.6, 0.1)
lambda_0 = lambda_0 * 1e-6
plt.plot(Omega_cut,smigma_cut)
plt.yscale('log')
plt.show()


pts = []
z = []

for i in lambda_0:
    etas = []
    for j,k in enumerate(sigma_cut):
        omega_0 = 2 * np.pi * c / i 
        x = HHG(2, 0.16*E_h/(e*a_0), omega_0, 100, 17, 1e30, k, 24.587387, 5*i/c)
        eta = x.eta_cutoff()
        etas.append(eta)
        pt = [i,Omega_cut[j],eta]
        pts.append(pt)
    z.append(etas)

#%%
lamb = []
sigma = []
eta = []
for i in pts:
    lamb.append(i[0])
    sigma.append(i[1])
    eta.append(i[2])

        
#%%
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
lamb, sigma = np.meshgrid(lambda_0, Omega_cut)
z = np.asarray(z)
# Plot the surface.
surf = ax.plot_surface(lamb, sigma, np.transpose(z), cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

#%%

I = 4e18
E_0 = np.sqrt(I * c * 4e-7 * np.pi)
lambda_0 = 8e-7
omega_0 = 2 * np.pi * c / lambda_0

t_rs_rs = []
t_rs_is = []
t_rl_rs = []
t_rl_is = []
t_is_rs = []
t_is_is = []
t_il_rs = []
t_il_is = []
harmonics = []
harmonic = np.arange(15,75.1,0.1)
for n in harmonic:
    try:
        x = HHG(10, E_0, omega_0, 100, n, 1e30, 4e-23, 21.5645, 5*lambda_0/c)
        t_rs, t_rl, t_is, t_il = x.saddle_solns_v1()
        t_rs_rs.append(t_rs.real)
        t_rs_is.append(t_rs.imag)
        t_rl_rs.append(t_rl.real)
        t_rl_is.append(t_rl.imag)
        t_is_rs.append(t_is.real)
        t_is_is.append(t_is.imag)
        t_il_rs.append(t_il.real)
        t_il_is.append(t_il.imag)
        harmonics.append(n)
    except:
        pass
#%%   
  

I = 4e18
E_0 = np.sqrt(I * c * 4e-7 * np.pi) / e
lambda_0 = 8e-7
omega_0 = 2 * np.pi * c / lambda_0

t_rs_rs = []
t_rs_is = []
t_rl_rs = []
t_rl_is = []
t_is_rs = []
t_is_is = []
t_il_rs = []
t_il_is = []
harmonics = []

for n in range(15,76,2):
    x = HHG(10, E_0, omega_0, 100, n, 1e30, 4e-23, 21.5645, 5*lambda_0/c)
    print(x._n_thresh)
    t_rs, t_rl, t_is, t_il = x.saddle_solns_v1()
    t_rs_rs.append(t_rs.real)
    t_rs_is.append(t_rs.imag)
    t_rl_rs.append(t_rl.real)
    t_rl_is.append(t_rl.imag)
    t_is_rs.append(t_is.real)
    t_is_is.append(t_is.imag)
    t_il_rs.append(t_il.real)
    t_il_is.append(t_il.imag)
    harmonics.append(n)

#%%
trsrs,trsis,harmonics2,tisrs,tisis = zip(*((x, y, z, a, b) for x, y, z, a, b\
                                    in zip(t_rs_rs,t_rs_is,harmonics,t_is_rs,t_is_is) \
                                                   if x <= 65 and x > 0 ))
trlrs,trlis,harmonics3,tilrs,tilis = zip(*((x,y,z,a,b) for x, y, z,a,b in \
                                     zip(t_rl_rs,t_rl_is,harmonics,t_il_rs,t_il_is)\
                               if x <= 88 and x > 62 ))
    
trsrs_arr = np.array(trsrs)
trsrs_arr = np.delete(trsrs_arr, -3)
trsis_arr = np.array(trsis)
trsis_arr = np.abs(trsis_arr)
trsis_arr = np.delete(trsis_arr, -3)
trsis_arr = -trsis_arr
harmonics2_arr = np.array(harmonics2)
harmonics2_arr = np.delete(harmonics2_arr, -3)
trlrs_arr = np.array(trlrs)
trlis_arr = np.array(trlis)
trlis_arr = np.abs(trlis_arr)
harmonics3_arr = np.array(harmonics3)

fig = plt.figure()
points = np.array([trsrs_arr, trsis_arr]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-2], points[1:-1], points[2:]], axis=1)

lc = LineCollection(segments, cmap=plt.get_cmap('Spectral'),
                        norm=plt.Normalize(15, 75))
lc.set_array(harmonics2_arr)
lc.set_linewidth(2)

points_l = np.array([trlrs_arr, trlis_arr]).T.reshape(-1, 1, 2)
segments_l = np.concatenate([points_l[:-2], points_l[1:-1], points_l[2:]], axis=1) 
lc2 = LineCollection(segments_l, cmap=plt.get_cmap('Spectral'),
                        norm=plt.Normalize(15, 75))
lc2.set_array(harmonics3_arr)
lc2.set_linewidth(2)

plt.gca().add_collection(lc)
plt.gca().add_collection(lc2)

axcb = fig.colorbar(lc)
axcb.set_label('Harmonic')
plt.xlim(0, 100)
plt.ylim(-20,20)
plt.xlabel('Real time of recombination (a.u.)')
plt.ylabel('Imaginary time of recombination (a.u.)')

plt.show()
#%%
plt.plot()
#%%
tisrs_arr = np.array(tisrs)
tisis_arr = np.array(tisis)
tisrs_arr = np.delete(tisrs_arr,-3)
tisis_arr = np.delete(tisis_arr,-3)
tilrs_arr = np.array(tilrs)
tilis_arr = np.array(tilis)
harmonics_is = np.array(harmonics2)
harmonics_is = np.delete(harmonics_is,-3)
harmonics_il = np.array(harmonics3)

fig = plt.figure()
points = np.array([tisrs_arr, tisis_arr]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-2], points[1:-1], points[2:]], axis=1)

lc = LineCollection(segments, cmap=plt.get_cmap('Spectral'),
                        norm=plt.Normalize(15, 75))
lc.set_array(harmonics_is)
lc.set_linewidth(2)

points_l = np.array([tilrs_arr, tilis_arr]).T.reshape(-1, 1, 2)
segments_l = np.concatenate([points_l[:-2], points_l[1:-1], points_l[2:]], axis=1) 
lc2 = LineCollection(segments_l, cmap=plt.get_cmap('Spectral'),
                        norm=plt.Normalize(15, 75))
lc2.set_array(harmonics_il)
lc2.set_linewidth(2)

plt.gca().add_collection(lc)
plt.gca().add_collection(lc2)

axcb = fig.colorbar(lc)
axcb.set_label('Harmonic')
plt.xlim(15, 30)
plt.ylim(-20,-5)
plt.xlabel('Real time of ionisation (a.u.)')
plt.ylabel('Imaginary time of ionisation (a.u.)')

plt.show()
#%%
E_0 = 0.125 * E_h / (e * a_0)
lambda_0 = 8.1e-7
omega_0 = 2 * np.pi * c / lambda_0
etas = []
harmonic = np.arange(40,77,0.1)
harmonics = []

for n in harmonic:
    try:
        x = HHG(2, E_0, omega_0, 100, n, 1e30, 4e-23, 24.587, 10)
        eta = x.eta_plateau()
        etas.append(eta)
        harmonics.append(n)
    except:
        pass
#%%
etas = np.array(etas)
#log_etas = np.log10(etas)
harmonics = np.array(harmonics)
#log_etas_fil = savgol_filter(log_etas,19,2)
#etas_fil = 10**(log_etas_fil)
etas_fil = savgol_filter(etas,11,2)
peaks, _ = find_peaks(etas_fil, height=0, prominence=1e-6, distance = 10 )
plt.plot(harmonics, etas_fil)
plt.yscale('log')
plt.xlabel('Harmonic')
plt.ylabel('Efficiency') 
plt.plot(harmonics[peaks], etas[peaks], "x")
plt.show()
#%%
E_0 = 0.151 * E_h / (e * a_0)
lambda_0 = 4.2e-7
omega_0 = 2 * np.pi * c / lambda_0
etas = []
harmonic = np.arange(13,20,0.01)
harmonics = []

for n in harmonic:
    try:
        x = HHG(2, E_0, omega_0, 100, n, 1e30, 4e-23, 24.587, 20*lambda_0/c)
        eta = x.eta_plateau()
        etas.append(eta)
        harmonics.append(n)
    except:
        pass

plt.plot(harmonics,etas)
plt.yscale('log')
plt.xlabel('Harmonic')
plt.ylabel('Efficiency') 
plt.show()  

#%%
E_0 = np.arange(0.1, 0.18, 0.0001)
lambda_0 = 8e-7
omega_0 = 2 * np.pi * c / lambda_0
etas = []
Es = []

for i in E_0:
    try:
        E = i * E_h / (e * a_0)
        x = HHG(10, E, omega_0, 0.1, 59, 1e30, 4e-23, 21.5645, 11)
        eta = x.eta_plateau()
        if eta < 1e-3:
            etas.append(eta)
            Es.append(i)
    except:
        pass
plt.plot(Es, etas)
plt.xlabel('E (a.u.)')
plt.ylabel('Efficiency')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.show()
#%%
etas_filtered = savgol_filter(etas,15,2)
plt.plot(Es, etas_filtered)
#%%
I = 4e18
E_0 = np.sqrt(I * c * 4e-7 * np.pi)
lambda_0 = 8e-7
omega_0 = 2 * np.pi * c / lambda_0
x = HHG(10,E_0, omega_0, 100, 17, 1e30,4e-23, 21.5645,10*lambda_0/c)
t = np.arange(0.1,50,0.1)
ass = []
for i in t:
    a = x.a_calc(i)
    ass.append(a)
#%%
plt.plot(t,ass)
plt.xlabel('time(a.u.)')
plt.ylabel('a(t)')
#%%
plt.subplot(2,2,1)
plt.plot(harmonics2_arr, trsrs_arr)
plt.xlabel('Harmonics')
plt.ylabel('Time (a.u.)')
plt.title('Recombination short')

plt.subplot(2,2,2)
plt.plot(harmonics3_arr, trlrs_arr)
plt.xlabel('Harmonics')
plt.ylabel('Time (a.u.)')
plt.title('Recombination long')

plt.subplot(2,2,3)
plt.plot(harmonics_is, tisrs_arr)
plt.xlabel('Harmonics')
plt.ylabel('Time (a.u.)')
plt.title('Ionisation short')

plt.subplot(2,2,4)
plt.plot(harmonics_il, tilrs_arr)
plt.xlabel('Harmonics')
plt.ylabel('Time (a.u.)')
plt.title('Ionisation long')

plt.tight_layout()
plt.show()

#%%
E_0_au = np.arange(0.01, 0.05, 0.0001)
lambda_0 = 8e-7
omega_0 = 2 * np.pi * c / lambda_0
etaps = []
etacs = []
Es = []
for i in E_0_au:
    try:
        E_0 = i * E_h / (e * a_0)
        x = HHG(10,E_0, omega_0, 100, 59, 1e30,4e-23, 21.5645,10*lambda_0/c)
        eta_p = x.eta_plateau()
        eta_c = x.eta_cutoff()
        etaps.append(eta_p)
        etacs.append(eta_c)
        Es.append(i)
    except:
        pass
plt.plot(Es, etaps)
plt.plot(Es, etacs)
plt.yscale('log')
plt.show()
#%%
I_p = 21.5645 * e / E_h
E_0 = 0.16 * E_h / (e * a_0)
lambda_0 = 8e-7
omega_0 = 2 * np.pi * c * hbar / (lambda_0 *  E_h)
harmonics = np.arange(15,76,0.1)
energies = []
a_rec_sq = []

def a_rec_calc(Z, k_n):
        a_0 = 1 / (2 * k_n * np.sqrt(np.pi))
        a_2 = -np.sqrt(5 / np.pi) / (2 * k_n)
        c_0 = np.sqrt(4 * np.pi / 3) * gaunt(1,1,0,0,0,0).n(10)
        c_2 = np.sqrt(4 * np.pi / 3) * gaunt(1,1,2,0,0,0).n(10)
        a_1 = 1j / (2 * k_n) * np.sqrt(3 / np.pi)
        c_1 = np.sqrt(4 * np.pi / 3) * gaunt(0,1,1,0,0,0).n(10)
        print(a_0, a_2, c_0, c_2, a_1, c_1)
        
        def u_01(r):
            return 2 / np.pi * k_n**2 * spherical_jn(0, k_n * r) * \
                        spherical_jn(1, k_n * r)
        def u_12(r):
            return 2 / np.pi * k_n**2 * spherical_jn(2, k_n * r) * \
                        spherical_jn(1, k_n * r)
        
        def u_01_r(r):
            return np.real(u_01(r))

        def u_01_i(r):
            return np.imag(u_01(r)) 
    
        def u_12_r(r):
            return np.real(u_12(r))
    
        def u_12_i(r):
            return np.imag(u_12(r))
        
        if Z == 2:
            if np.iscomplex(k_n):
                I_01_r = integrate.quad(u_01_r,0,np.inf)
                I_01_i = integrate.quad(u_01_i,0,np.inf)
                I_01 = I_01_r[0] + 1j * I_01_i[0]
            else:
                I_01 = integrate.quad(u_01, 0, np.inf)[0]
            print(I_01)
            return -a_1 * c_1 * 2 * I_01
        else:
            if np.iscomplex(k_n):
                I_01_r = integrate.quad(u_01_r,0,np.inf)
                I_01_i = integrate.quad(u_01_i,0,np.inf)
                I_12_r = integrate.quad(u_12_r,0,np.inf)
                I_12_i = integrate.quad(u_12_i,0,np.inf)
                I_01 = I_01_r[0] + 1j * I_01_i[0]
                I_12 = I_12_r[0] + 1j * I_12_i[0]
            else:
                I_01 = integrate.quad(u_01, 0, np.inf)[0]
                I_12 = integrate.quad(u_12, 0, np.inf)[0]
            return -a_0 * c_0 * Z * I_01 - a_2 * c_2 * Z * I_12

for n in harmonics:
    E = n * omega_0 * E_h / e
    energies.append(E)
    k_n = cmath.sqrt(2 * (n * omega_0 - I_p))
    a = a_rec_calc(10,k_n)
    a = abs(a)**2
    a_rec_sq.append(a)

plt.plot(energies, a_rec_sq) 
plt.yscale('log')  

#%%
default_params = {
    'Atomic Number' : 10,
    'Laser E-field Amplitude' : 7.2e10,
    'Laser Angular Frequency' : 2354564459136066.5,
    'Max Harmonic' : 75,
    'Laser Cycles' : 5,
    'Propagation Distance' : 100,
    'Number Density' : 5e25,
    'Abs. Cross-section' : 4e-23
    }
sample = json.dumps(sample_params, indent=4)

with open("sample_params.json", "w") as outfile:
    outfile.write(sample)
#%%
def writeToJSONFile(path, fileName, data):
        json.dump(data, path)

window = Tk()
window.geometry('500x640')
window.title('Input')

def check():
    params = {
        'Atomic Number' : 10,
        'Laser E-field Amplitude' : 7.2e10,
        'Laser Angular Frequency' : 2354564459136066.5,
        'Max Harmonic' : 75,
        'Laser Cycles' : 5,
        'Propagation Distance' : 100,
        'Number Density' : 5e25,
        'Abs. Cross-section' : 4e-23
    }
    params_list = list(params)
    a = Z.get()
    b = E_0.get()
    c = Omega_0.get()
    d = Max_n.get()
    e = Cycles.get()
    f = L.get()
    g = N_density.get()
    h = Sigma.get()
    entries = [a,b,c,d,e,f,g,h]
    for i,j in enumerate(entries):
        if type(j) == int or type(j) == float:
            params[params_list[i]] = j
        else:
            pass
    files = [('JSON File', '*.json')]
    fileName='params'
    filepos = asksaveasfile(filetypes = files,defaultextension = json,initialfile='params')
    writeToJSONFile(filepos, fileName, params)

z = Label(window, text='Atomic Number:')
Z = Entry(window)
e_0 = Label(window, text = 'Laser Electric Field Amplitude (V/m):')
E_0 = Entry(window)
omega_0 = Label(window, text = 'Laser Angular Frequency (rad/sec):')
Omega_0 = Entry(window)
max_n = Label(window, text = 'Max Harmonic:')
Max_n = Entry(window)
cycles = Label(window, text = 'Laser Cycles:')
Cycles = Entry(window)
l = Label(window, text = 'Propagation Distance (m):')
L = Entry(window)
n_density = Label(window, text = 'Number Density (m^-3):')
N_density = Entry(window)
sigma = Label(window, text = 'Absorption Cross-section (m^2):')
Sigma = Entry(window)
submit = Button(window,text='Submit',command = check).grid(row=8, column=1)

z.grid(row=0, column=0)
e_0.grid(row=1, column=0)
omega_0.grid(row=2, column=0)
max_n.grid(row=3, column=0)
cycles.grid(row=4, column=0)
l.grid(row=5, column=0)
n_density.grid(row=6, column=0)
sigma.grid(row=7, column=0)
Z.grid(row=0, column=1)
E_0.grid(row=1, column=1)
Omega_0.grid(row=2, column=1)
Max_n.grid(row=3, column=1)
Cycles.grid(row=4, column=1)
L.grid(row=5, column=1)
N_density.grid(row=6, column=1)
Sigma.grid(row=7, column=1)


mainloop()

