# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 17:49:04 2021

@author: joshu
"""
import numpy as np
import cmath
import math
import scipy.integrate as integrate
import mpmath
from sympy.physics.wigner import gaunt
from scipy.special import spherical_jn
from scipy.optimize import newton
import types

#constants
alpha = 1/137
a_0 = 5.29177210903e-11
E_h = 4.3597447222071e-18
c = 299792458
hbar =  6.62607015e-34/(2 * np.pi)
e = 1.602176634e-19
m_e = 9.1093837015e-31

def gamma_integrand(x,a): #gamma function for calculation of ionisation rate, w
    return x**(a-1)*np.exp(-x)

def gamma_func_calc(z): 
    return integrate.quad(gamma_integrand,0,np.inf,args=(z))[0]

class HHG:
    def __init__(self, Z, E_0, omega_0, L, harmonic, n_density, sigma, I_p, cycles):
        #defining variables and converting to atomic units
        #E_0, I_p in eV
        self._Z = Z
        self._E_0 = E_0 * e * a_0 / E_h  
        self._omega_0 = omega_0 * hbar / E_h
        self._harmonic = harmonic
        self._n_density = n_density * a_0**3
        self._omega_p = 2 * np.sqrt(np.pi * n_density)
        self._sigma = sigma / a_0**2
        self._I_p = I_p * e / E_h
        self._L = L / a_0
        self._L_abs = 1/(n_density * sigma)
        self._delta_k = -harmonic * self._omega_0 * self._omega_p**2 * alpha / \
            (2 * c * self._omega_0**2) 
        self._T = 2 * np.pi * cycles * E_h / (hbar * omega_0) 
        self._N = cycles
        self._Omega = harmonic * self._omega_0
        self._U_p = (self._E_0 / (2 * self._omega_0))**2
        self._Omega_cut = self._I_p + 3.17 * self._U_p
        self._tb_cut = 1.884 / self._omega_0
        self._ta_cut = 5.97 / self._omega_0
        self._sigma_cut = self._sigma
        self._beta = self.beta_calc()
        self._kappa_0 = abs(self.a_calc(self._tb_cut) * \
                            self.a_calc(self._ta_cut))**2
        self._p_z = alpha * self._omega_0
        self._Keldysh = np.sqrt(0.5 * self._I_p / self._U_p)
        self._k_n = cmath.sqrt(2 * (self._harmonic * self._omega_0 - self._I_p))
        self._a_rec = self.a_rec_calc()
        if L > 3 * self._L_abs: #prop. distance above which g can be approximated as 1
            self._g = 1
        else: 
            self._g = cmath.exp(self._delta_k * self._L * 1j) - np.exp(-0.5 * self._L / \
            self._L_abs) / (1 + 2j*(self._delta_k * self._L_abs))
        self._v_0 = self._E_0 / self._omega_0
        self._Gamma_n = np.sqrt(2 * (self._harmonic * self._omega_0 - self._I_p))
        self._Gamma = -np.sqrt(2 * self._I_p) * 1j
        self._n_thresh = math.ceil(self._I_p / self._omega_0)
        if harmonic < self._n_thresh:
            print('Chosen harmonic is below the threshold for ionisation')
        self._t_rs, self._t_rl, self._t_is, self._t_il = self.saddle_solns_v2()
        self._t_rs = self._t_rs.real
        self._t_rl = self._t_rl.real
        self._t_is = self._t_is.real
        self._t_il = self._t_il.real
        
        
    def E_driving(self, t): #driving field of laser
        if t < 0 or t > self._T:
            return 0
        else:
            return self._E_0 * np.sin(self._omega_0 * t)
        
    def w_calc(self,t):
        E = abs(self.E_driving(t))
        n = 1 / np.sqrt(2 * self._I_p)
        C_nl_sq = (2 * np.exp(1) / n)**(2 * n) / (2 * np.pi * n)
        if self._Z == 2:
            f = 1
        else:
            f = 3
        return C_nl_sq * f * self._I_p * (2 * np.sqrt((2 * self._I_p)**3) / abs(E))**(2 * n - 1) \
                * np.exp(- 2 * np.sqrt((2 * self._I_p)**3) / (3 * abs(E)))               
  
    def beta_calc(self):
        I = integrate.quad(self.w_calc, 0, np.pi / self._omega_0)
        return np.exp(-I[0])
    
    def a_calc(self, tt): 
        I = integrate.quad(self.w_calc, 0, tt)
        return np.exp(-0.5 * I[0])
        
    def eta_cutoff(self):
        f = np.sqrt(2 * self._I_p) * self._omega_0**5 * self._a_rec**2 \
            * abs(self._g)**2 / (self._E_0**(16 / 3) * self._Omega_cut**2 * \
            self._sigma_cut**2)
        b = (1 - self._beta**(4 * (self._N - 1)) * abs(1 + self._beta)**2) / \
            ((1 - self._beta**4) * self._N)
        w = self.w_calc(self._tb_cut)
        return 0.0236 * f * b * w * self._kappa_0
            
    def A_potential(self, t):
        if isinstance(t, complex):
            return self._E_0 * cmath.cos(self._omega_0 * t) / self._omega_0
        else:
            return self._E_0 * np.cos(self._omega_0 * t) / self._omega_0
    
    def eta_plateau(self):
        f1 = np.sqrt(2 * self._I_p) * self._omega_0**5 * abs(self._a_rec)**2 \
             * abs(self._g)**2 / (self._E_0**4 * self._Omega**2 * self._sigma**2)
        b = (1 - self._beta**(4 * (self._N - 1))) / ((1 - self._beta**4) * self._N) \
            * abs(1 + self._beta * cmath.exp(np.pi * (1 - self._Omega / self._omega_0) \
            * 1j))**2
        fs = self.a_calc(self._t_is) * self.a_calc(self._t_rs) * np.sqrt( \
             self.w_calc(self._t_is)) / (cmath.sin(self._omega_0 * \
             self._t_is) * (self._omega_0 * (self._t_rs - self._t_is) / (2 * np.pi))**1.5)\
             * cmath.exp(-1j * (self.S_action(self._t_is, self._t_rs) - self._Omega\
             * self._t_rs)) / np.sqrt(abs(self.S_2nd_der(self._t_is, self._t_rs))) 
        fl = self.a_calc(self._t_il) * self.a_calc(self._t_rl) * np.sqrt( \
             self.w_calc(self._t_il)) / (cmath.sin(self._omega_0 * \
             self._t_il) * (self._omega_0 * (self._t_rl - self._t_il) / (2 * np.pi))**1.5)\
             * cmath.exp(-1j * (self.S_action(self._t_il, self._t_rl) - self._Omega\
             * self._t_rl - np.pi / 2)) / \
             np.sqrt(abs(self.S_2nd_der(self._t_il, self._t_rl)))
        print(f1,b,fs,fl)
        return 0.0107 * f1 * b * abs(fs + fl)**2
    
    def a_rec_calc(self):
        a_0 = 1 / (2 * self._k_n * np.sqrt(np.pi))
        a_2 = -np.sqrt(5 / np.pi) / (2 * self._k_n)
        c_0 = np.sqrt(4 * np.pi / 3) * gaunt(1,1,0,0,0,0).n(10)
        c_2 = np.sqrt(4 * np.pi / 3) * gaunt(1,1,2,0,0,0).n(10)
        a_1 = 1j / (2 * self._k_n) * np.sqrt(3 / np.pi)
        c_1 = np.sqrt(4 * np.pi / 3) * gaunt(0,1,1,0,0,0).n(10)
        print(a_0, a_2, c_0, c_2, a_1, c_1)
        
        def u_01(r):
            return 2 / np.pi * self._k_n**2 * spherical_jn(0, self._k_n * r) * \
                        spherical_jn(1, self._k_n * r)
        def u_12(r):
            return 2 / np.pi * self._k_n**2 * spherical_jn(2, self._k_n * r) * \
                        spherical_jn(1, self._k_n * r)
        
        def u_01_r(r):
            return np.real(u_01(r))

        def u_01_i(r):
            return np.imag(u_01(r)) 
    
        def u_12_r(r):
            return np.real(u_12(r))
    
        def u_12_i(r):
            return np.imag(u_12(r))
        
        if self._Z == 2:
            if np.iscomplex(self._k_n):
                I_01_r = integrate.quad(u_01_r,0,np.inf)
                I_01_i = integrate.quad(u_01_i,0,np.inf)
                I_01 = I_01_r[0] + 1j * I_01_i[0]
            else:
                I_01 = integrate.quad(u_01, 0, np.inf)[0]
            print(I_01)
            return -a_1 * c_1 * 2 * I_01
        else:
            if np.iscomplex(self._k_n):
                I_01_r = integrate.quad(u_01_r,0,np.inf)
                I_01_i = integrate.quad(u_01_i,0,np.inf)
                I_12_r = integrate.quad(u_12_r,0,np.inf)
                I_12_i = integrate.quad(u_12_i,0,np.inf)
                I_01 = I_01_r[0] + 1j * I_01_i[0]
                I_12 = I_12_r[0] + 1j * I_12_i[0]
            else:
                I_01 = integrate.quad(u_01, 0, np.inf)[0]
                I_12 = integrate.quad(u_12, 0, np.inf)[0]
            return -a_0 * c_0 * self._Z * I_01 - a_2 * c_2 * self._Z * I_12
        
    def p_mom(self, t_i, t, v_i=0):
        return v_i - self.A_potential(t_i) + 2 * self.A_potential(t) 
    
    def S_integrand(self, t_i, t):
        return (self.p_mom(t_i, t) - self.A_potential(t_i))**2
        
    def S_action(self, t1, t2):
        I = integrate.quad(self.S_integrand, t1, t2, args=(t2))[0]
        return 0.5 * I + self._I_p * (t2 - t1)
    
    def S_2nd_der(self, t_i, t_r):
        return (self.A_potential(t_i) - self.p_mom(t_i, t_r)) * self.E_driving(t_i)
    
    def t_i(self, t_r):
        p_r = self._Gamma_n + self.A_potential(t_r)
        return 1 / self._omega_0 * cmath.acos(1 / self._v_0 * (p_r - self._Gamma))
    
    def saddle_eq_tr(self, t_r):
        t_i = self.t_i(t_r)
        return (self._Gamma_n + self.A_potential(t_r)) * (t_r - t_i) + self._v_0 / \
            self._omega_0 * (cmath.sin(self._omega_0 * t_r) - \
                             cmath.sin(self._omega_0 * t_i))
    def saddle_eq_tr_prime(self, t_r):
        p_r = self._Gamma_n + self.A_potential(t_r)
        t_i = self.t_i(t_r)
        dt_i = cmath.sin(self._omega_0 * t_r) / cmath.sqrt(1 - ((p_r - self._Gamma)\
                                                                / self._v_0)**2)
        return (t_r - t_i) * (-self._E_0 * cmath.sin(self._omega_0 * t_r)) - \
               (self._Gamma_n + self.A_potential(t_r)) * dt_i + 2 * \
               self.A_potential(t_r) - self.A_potential(t_i) * dt_i + self._Gamma_n
    
    def set_harmonic(self, new_n):
        self._harmonic = new_n
    
    def saddle_solns_v2(self):
        harmonic = self._harmonic
        self.set_harmonic(self._n_thresh)
        T_0 = 2 * np.pi / self._omega_0
        guess_1 = 0.025 * T_0 + 0.1j * T_0
        t_rs_1 = newton(self.saddle_eq_tr, guess_1, fprime=self.saddle_eq_tr_prime,\
                        maxiter=200)
        if t_rs_1.imag > 0:
            t_rs_1 = np.conj(t_rs_1)
        guess_2 = t_rs_1 + T_0 / 2
        t_rl_1 = newton(self.saddle_eq_tr, guess_2, fprime=self.saddle_eq_tr_prime,\
                        maxiter=200)
        print(t_rs_1, t_rl_1)
        self.set_harmonic(harmonic)
        t_rs = newton(self.saddle_eq_tr, t_rs_1, fprime=self.saddle_eq_tr_prime)
        if t_rs.imag > 0:
            t_rs = np.conj(t_rs)
        t_rl = newton(self.saddle_eq_tr, t_rl_1, fprime=self.saddle_eq_tr_prime)
        t_is = self.t_i(t_rs)
        t_il = self.t_i(t_rl)
        return t_rs, t_rl, t_is, t_il
    
    def saddle_solns_v1(self):
        T_0 = 2 * np.pi / self._omega_0
        guess_1 = 0.025 * T_0 + 0.1j * T_0
        t_rs = newton(self.saddle_eq_tr, guess_1, fprime=self.saddle_eq_tr_prime)
        #if t_rs.imag > 0:
         #   t_rs = np.conj(t_rs)
        guess_2 = t_rs + T_0 / 2
        t_rl = newton(self.saddle_eq_tr, guess_2, fprime=self.saddle_eq_tr_prime)
        t_is = self.t_i(t_rs)
        t_il = self.t_i(t_rl)
        return t_rs, t_rl, t_is, t_il
    
    
        
        
        
        
        
        
        
        