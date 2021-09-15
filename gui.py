# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:35:49 2021

@author: joshu
"""
from tkinter import *
import json
from tkinter.filedialog import asksaveasfile
import glob
import os

def writeToJSONFile(path, data):
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
        try:
            value = float(j)
            params[params_list[i]] = value
        except:
            pass
    files = [('JSON File', '*.json')]
    filepos = asksaveasfile(filetypes = files,defaultextension = json,initialfile='params')
    writeToJSONFile(filepos, params)
    window.destroy()


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

LatestFile = max(glob.iglob('*.json'),key=os.path.getmtime)
