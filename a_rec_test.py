# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 11:07:43 2021

@author: joshu
"""
from hhg import *
import matplotlib.pyplot as plt

I_p = I_ps[10] * e / E_h
E_0 = 0.16 * E_h / (e * a_0)
lambda_0 = 8e-7
omega_0 = 2 * np.pi * c * hbar / (lambda_0 *  E_h)
harmonics = np.arange(1,200,0.1)
energies = []
a_rec_sq = []
a_s = []
I_01s = []
I_12s = []

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
                I_01 = integrate.quad(u_01, 0, 100)[0]
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
                I_01 = integrate.quad(u_01, 0, 100)[0]
                I_12 = integrate.quad(u_12, 0, 100)[0]
            return -a_0 * c_0 * Z * I_01 - a_2 * c_2 * Z * I_12, I_01, I_12

for n in harmonics:
    E = n * omega_0 * E_h / e
    energies.append(E)
    k_n = cmath.sqrt(2 * (n * omega_0 - I_p))
    a, I_01, I_12 = a_rec_calc(10,k_n)
    a_sq = abs(a)**2
    a_rec_sq.append(a_sq)
    a_s.append(abs(a))
    I_01s.append(I_01)
    I_12s.append(I_12)
#%%
fig, axs = plt.subplots(2,2)
axs[0,0].plot(energies, a_s) 
axs[0,0].yscale('log')
axs[0,0].xlabel('Energy (eV)')
axs[0,0].ylabel('|a_rec|')
axs[0,1].plot(energies, I_01s)
axs[0,1].yscale('log')
axs[1,0].plot(energies, I_12s)
axs[1,0].yscale('log')
plt.show()
#%%
k_n = cmath.sqrt(2 * (14*omega_0 - I_p))
def u_00(r):
    return 2 / np.pi * k_n**2 * spherical_jn(0, k_n * r) * \
                        spherical_jn(1, k_n * r)
def u_02(r):
    return 2 / np.pi * k_n**2 * spherical_jn(2, k_n * r) * \
                        spherical_jn(1, k_n * r)
r = np.arange(0,100,0.1)

u00 = [u_00(x) for x in r]
u02 = [u_02(x) for x in r]
plt.plot(r,u00)
plt.plot(r,u02)

#%%

I_p = 24.587387 * e / E_h
E_0 = 0.16 * E_h / (e * a_0)
lambda_0 = 8e-7
omega_0 = 2 * np.pi * c * hbar / (lambda_0 *  E_h)
harmonics = np.arange(1,100,0.1)
energies = []
a_rec_sq = []
a_s = []

def a_rec_calc(Z, k_n):
        a_0 = 1 / (2 * k_n * np.sqrt(np.pi))
        a_2 = -np.sqrt(5 / np.pi) / (2 * k_n)
        c_0 = np.sqrt(4 * np.pi / 3) * gaunt(1,1,0,0,0,0).n(10)
        c_2 = np.sqrt(4 * np.pi / 3) * gaunt(1,1,2,0,0,0).n(10)
        a_1 = 1j / (2 * k_n) * np.sqrt(3 / np.pi)
        c_1 = np.sqrt(4 * np.pi / 3) * gaunt(0,1,1,0,0,0).n(10)
        print(a_0, a_2, c_0, c_2, a_1, c_1)
        
        def u_01(r):
            return 2 / np.pi * k_n**2 * r**3 * spherical_jn(0, k_n * r) * \
                        spherical_jn(1, k_n * r)
        def u_12(r):
            return 2 / np.pi * k_n**2 * spherical_jn(2, k_n * r) * \
                        spherical_jn(0, k_n * r)
        
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
                I_01 = integrate.quad(u_01, 0, 100)[0]
            print(I_01)
            return a_1 * c_1 * I_01
        else:
            if np.iscomplex(k_n):
                I_01_r = integrate.quad(u_01_r,0,np.inf)
                I_01_i = integrate.quad(u_01_i,0,np.inf)
                I_12_r = integrate.quad(u_12_r,0,np.inf)
                I_12_i = integrate.quad(u_12_i,0,np.inf)
                I_01 = I_01_r[0] + 1j * I_01_i[0]
                I_12 = I_12_r[0] + 1j * I_12_i[0]
            else:
                I_01 = integrate.quad(u_01, 0, 100)[0]
                I_12 = integrate.quad(u_12, 0, 100)[0]
            return -a_0 * c_0 * Z * I_01 - a_2 * c_2 * Z * I_12

for n in harmonics:
    E = n * omega_0 * E_h / e
    energies.append(E)
    k_n = cmath.sqrt(2 * (n * omega_0 - I_p))
    a = a_rec_calc(2,k_n) 
    a_sq = abs(a)**2
    a_rec_sq.append(a_sq)
    a_s.append(abs(a))

plt.plot(energies, a_s) 
plt.yscale('log') 
