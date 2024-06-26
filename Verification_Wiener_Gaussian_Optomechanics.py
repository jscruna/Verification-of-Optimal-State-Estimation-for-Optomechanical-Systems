#THIS IS THE FINAL VERSION OF THE CODE

#function to integrate over an interval
#It returns the frequency in the interval and the sliced function
def bandpower(pxx,freq,freq_int):
    import numpy as np
    # Indices corresponding to 10 Hz and 50 Hz
    start_index = np.searchsorted(freq, freq_int[0])
    end_index = np.searchsorted(freq, freq_int[1])
    # Slicing the arrays
    frequency_slice = freq[start_index:end_index]
    psd_slice = pxx[start_index:end_index]
    #psd_slice = Sqq(2*np.pi*frequency_slice)
    # Integrating over the specified range
    return 2*np.trapz(psd_slice, 2*np.pi*frequency_slice)/(2*np.pi)

def pxx_fun(data,fs,fmin):
    import numpy as np
    from scipy.signal import welch
    import math
    """Calculate the nearest integer of the ratio between frequency sample and frequency resolution"""
    ratio= math.ceil(fs/fmin)# this ratio determines the vector size 
    #nperseg_val = pow2_nextpow2(ratio)
    nperseg_val = 2 ** (math.ceil(math.log(ratio, 2)))
    f, Pxx_den = welch(data, fs, window='hann', nperseg=nperseg_val, noverlap=None)
    return f, Pxx_den

import numpy as np
import scipy.constants as const
from scipy.signal import butter, filtfilt
from scipy.integrate import trapz
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import h5py
import pandas as pd
import numpy as np
# Number of samples for the Monte Carlo simulation
sam = 100

# Generate different seeds to avoid dependence between random variables
np.random.seed(1)
rand1 = np.random.randn(sam)
np.random.seed(2)
rand2 = np.random.randn(sam)
np.random.seed(3)
rand3 = np.random.randn(sam)
np.random.seed(4)
rand4 = np.random.randn(sam)
np.random.seed(5)
rand5 = np.random.randn(sam)
np.random.seed(6)
rand6 = np.random.randn(sam)
np.random.seed(7)
rand7 = np.random.randn(sam)

# Optomechanical parameters
mass = 7.71e-6 + 0.01 * rand1 * 1e-6
mass_mean = 7.71e-6

kappa_p = 2 * np.pi * (1.64 + 0.02 * rand2) * 1e6
kappa_p_mean = 2 * np.pi * 1.64 * 1e6

g = -2 * np.pi * (3.2 + 0.2 * rand3) * 1e4
g_mean = -2 * np.pi * 3.2 * 1e4

omega_m = 2 * np.pi * (280 + 7 * rand4)
omega_m_mean = 2 * np.pi * 280

Detune = 0.0292 + 0.0004 * rand5
Detune_mean = 0.0292

Temp = (11 + 2 * rand6) * 1e-3
Temp_mean = 11e-3

n_th = (8 + 2 * rand7) * 1e5
n_th_mean = 8e5

gamma_m = 1.1 * 2 * np.pi
gamma_m_mean = 1.1 * 2 * np.pi

finesse = 1849
beta = 0.631
N_th = 19
hbar = const.hbar
k_B = const.Boltzmann
c = const.c
lambda_ = 1064e-9
omega_laser = c / lambda_ * 2 * np.pi
P_in = 30e-3 * 0.8453
L = 0.098
eta = 0.92
kappa_in_over_all = 0.2247

x_zpf = np.sqrt(hbar / 2 / mass / omega_m)
x_zpf_mean = np.sqrt(hbar / 2 / mass_mean / omega_m_mean)
p_zpf_mean = np.sqrt(mass_mean * omega_m_mean * hbar / 2)

# Wiener Filter parameters
n_x = 2 * gamma_m * (2 * n_th + 1) + 16 * g ** 2 * (2 * N_th + 1) / ((1 + 4 * Detune ** 2) * kappa_p)
n_x_mean = 2 * gamma_m * (2 * n_th_mean + 1) + 16 * g_mean ** 2 * (2 * N_th + 1) / ((1 + 4 * Detune_mean ** 2) * kappa_p_mean)
lambda_x = 64 * g ** 2 * eta * Detune ** 2 / ((2 * eta * N_th + 1) * kappa_p * (1 + 4 * Detune ** 2) ** 2)
lambda_x_mean = 64 * g_mean ** 2 * eta * Detune_mean ** 2 / ((2 * eta * N_th + 1) * kappa_p_mean * (1 + 4 * Detune_mean ** 2) ** 2)
sigma_x = -32 * g ** 2 * eta * Detune * (2 * N_th + 1) / ((1 + 4 * Detune ** 2) ** 2 * kappa_p * (2 * eta * N_th + 1))
sigma_x_mean = -32 * g_mean ** 2 * eta * Detune_mean * (2 * N_th + 1) / ((1 + 4 * Detune_mean ** 2) ** 2 * kappa_p_mean * (2 * eta * N_th + 1))
omega_x = np.sqrt(np.sqrt(omega_m ** 4 + 2 * sigma_x * omega_m ** 3 + n_x * lambda_x * omega_m ** 2))
omega_x_mean = np.sqrt(np.sqrt(omega_m_mean ** 4 + 2 * sigma_x_mean * omega_m_mean ** 3 + n_x_mean * lambda_x_mean * omega_m_mean ** 2))
gamma_x = np.sqrt(gamma_m ** 2 - 2 * omega_m * (omega_m + sigma_x) + 2 * omega_x ** 2)
gamma_x_mean = np.sqrt(gamma_m ** 2 - 2 * omega_m_mean * (omega_m_mean + sigma_x_mean) + 2 * omega_x_mean ** 2)
M = (2 * N_th * eta + 1)

# Theoretical values of the variances for prediction and retrodiction
V11_aux = (gamma_x_mean - gamma_m) / lambda_x_mean
VE11_aux = (gamma_x_mean + gamma_m) / lambda_x_mean

V22_aux = (gamma_x_mean - gamma_m) / (lambda_x_mean * 2 * omega_m_mean ** 2) * (2 * omega_m_mean * (omega_m_mean + sigma_x_mean) + gamma_x_mean * (gamma_x_mean - gamma_m))
VE22_aux = (gamma_x_mean + gamma_m) / (lambda_x_mean * 2 * omega_m_mean ** 2) * (2 * omega_m_mean * (omega_m_mean + sigma_x_mean) + gamma_x_mean * (gamma_x_mean + gamma_m))

V12_aux = (gamma_x_mean - gamma_m) ** 2 / (2 * omega_m_mean * lambda_x_mean)
VE12_aux = -(gamma_x_mean + gamma_m) ** 2 / (2 * omega_m_mean * lambda_x_mean)

# Verification process under different prediction and retrodiction filters
resol = np.linspace(1, 120, 120)
var_test = np.zeros((len(resol), sam))
var_test_p = np.zeros((len(resol), sam))
var_test_qp = np.zeros((len(resol), sam))

#THIS IS FOR TESTING THE NEW FUNCTION TO OBTAIN A BETTER SAMPLE
def process_sample(j):
    fmin = 1
    rPD2 = pd.read_csv('your_experimental_data.csv',header=None,names=['V'])
    rPD = rPD2.values.flatten()
    for n_notch in range(1, 41):
        b, a = signal.butter(1, [(n_notch + 1) * 50 - fmin/2, (n_notch + 1) * 50 + fmin/2], btype='stop', fs=fs)
        rPD = signal.filtfilt(b, a, rPD)
    rPD3 = rPD
    cal = -np.pi * c * mass[j] / finesse / np.cos(beta) * (1 - kappa_in_over_all) * (omega_m[j] ** 2 * 0.73 * 500)
    rPD3 /= cal
    rPD3 /= x_zpf[j]
    obs2 = rPD3 * np.sqrt(lambda_x[j] * M)

    Hq_pred = lambda omega: (1 / np.sqrt(lambda_x[j] * M)) * ((omega_x[j] ** 2 - omega_m[j] ** 2) - 1j * omega * (gamma_x[j] - gamma_m)) / ((omega_x[j] ** 2 - omega ** 2) + 1j * gamma_x[j] * omega)
    Hq_retr = lambda omega: (1 / np.sqrt(lambda_x[j] * M)) * ((omega_x[j] ** 2 - omega_m[j] ** 2) + 1j * omega * (gamma_x[j] + gamma_m)) / np.conj((omega_x[j] ** 2 - omega ** 2) + 1j * gamma_x[j] * omega)
    Hp_pred = lambda omega: 1 / (np.sqrt(lambda_x[j] * M) * omega_m[j]) * (-(gamma_x[j] - gamma_m) * omega_m[j] ** 2 - 1j * omega * (omega_x[j] ** 2 - omega_m[j] ** 2 + (gamma_x[j] - gamma_m) ** 2 - (gamma_x[j] - gamma_m) * gamma_x[j])) / ((omega_x[j] ** 2 - omega ** 2) + 1j * gamma_x[j] * omega)
    Hp_retr = lambda omega: 1 / (np.sqrt(lambda_x[j] * M) * omega_m[j]) * ((gamma_x[j] + gamma_m) * omega_m[j] ** 2 - 1j * omega * (omega_x[j] ** 2 - omega_m[j] ** 2 + (gamma_x[j] + gamma_m) ** 2 - (gamma_x[j] + gamma_m) * gamma_x[j])) / np.conj((omega_x[j] ** 2 - omega ** 2) + 1j * gamma_x[j] * omega)

    var_test_val = np.zeros(len(resol))
    var_test_p_val = np.zeros(len(resol))
    var_test_qp_val = np.zeros(len(resol))
    #results = np.zeros((len(resol), 3))
    for mm in range(len(resol)):
        fmin = resol[mm]
        Freq5, Pxx5 = pxx_fun(obs2, fs, fmin)
        #var_test_val = trapz(Hq_retr(2 * np.pi * Freq5) * np.conj(Hq_retr(2 * np.pi * Freq5)) * Pxx5 - Hq_pred(2 * np.pi * Freq5) * np.conj(Hq_pred(2 * np.pi * Freq5)) * Pxx5, Freq5)
        #var_test_p_val = trapz(Hp_retr(2 * np.pi * Freq5) * np.conj(Hp_retr(2 * np.pi * Freq5)) * Pxx5 - Hp_pred(2 * np.pi * Freq5) * np.conj(Hp_pred(2 * np.pi * Freq5)) * Pxx5, Freq5)
        #var_test_qp_val = trapz(np.real(Hq_retr(2 * np.pi * Freq5) * np.conj(Hp_retr(2 * np.pi * Freq5)) * Pxx5 - Hq_pred(2 * np.pi * Freq5) * np.conj(Hp_pred(2 * np.pi * Freq5)) * Pxx5), Freq5)
        #var_test_val = np.real(bandpower(Hq_retr(2 * np.pi * Freq5) * np.conj(Hq_retr(2 * np.pi * Freq5)) * Pxx5 - Hq_pred(2 * np.pi * Freq5) * np.conj(Hq_pred(2 * np.pi * Freq5)) * Pxx5,Freq5,[200,1000]))
        #var_test_p_val = np.real(bandpower(Hp_retr(2 * np.pi * Freq5) * np.conj(Hp_retr(2 * np.pi * Freq5)) * Pxx5 - Hp_pred(2 * np.pi * Freq5) * np.conj(Hp_pred(2 * np.pi * Freq5)) * Pxx5,Freq5,[200,1000]))
        #var_test_qp_val = np.real(bandpower(Hq_retr(2 * np.pi * Freq5) * np.conj(Hp_retr(2 * np.pi * Freq5)) * Pxx5 - Hq_pred(2 * np.pi * Freq5) * np.conj(Hp_pred(2 * np.pi * Freq5)) * Pxx5,Freq5,[200,1000]))
        var_test_val[mm] = np.real(bandpower(Hq_retr(2 * np.pi * Freq5) * np.conj(Hq_retr(2 * np.pi * Freq5)) * Pxx5 - Hq_pred(2 * np.pi * Freq5) * np.conj(Hq_pred(2 * np.pi * Freq5)) * Pxx5,Freq5,[150,1000]))/2
        var_test_p_val[mm] = np.real(bandpower(Hp_retr(2 * np.pi * Freq5) * np.conj(Hp_retr(2 * np.pi * Freq5)) * Pxx5 - Hp_pred(2 * np.pi * Freq5) * np.conj(Hp_pred(2 * np.pi * Freq5)) * Pxx5,Freq5,[150,1000]))/2
        var_test_qp_val[mm] = np.real(bandpower(Hq_retr(2 * np.pi * Freq5) * np.conj(Hp_retr(2 * np.pi * Freq5)) * Pxx5 - Hq_pred(2 * np.pi * Freq5) * np.conj(Hp_pred(2 * np.pi * Freq5)) * Pxx5,Freq5,[150,1000]))
        #results[mm, :] = [var_test_val, var_test_p_val, var_test_qp_val]

    return var_test_val,var_test_p_val

#process_sample(1)
results = Parallel(n_jobs=5)(delayed(process_sample)(j) for j in range(sam))
#results = np.array(results)
import numpy as np

#results = np.array(results)

x = resol
y_mean = np.zeros(len(resol))
y_p_mean = np.zeros(len(resol))
#y_qp_mean = np.zeros(len(resol))
y_upper = np.zeros(len(resol))
y_lower = np.zeros(len(resol))
y_upper_p = np.zeros(len(resol))
y_lower_p = np.zeros(len(resol))
#y_upper_qp = np.zeros(len(resol))
#y_lower_qp = np.zeros(len(resol))

for kk in range(len(resol)):
    y_upper[kk] = np.mean(np.array(results).transpose()[kk,0]) + np.std(np.array(results).transpose()[kk,0])
    y_lower[kk] = np.mean(np.array(results).transpose()[kk,0]) - np.std(np.array(results).transpose()[kk,0])
    y_mean[kk] = np.mean(np.array(results).transpose()[kk,0])
    y_upper_p[kk] = np.mean(np.array(results).transpose()[kk,1]) + np.std(np.array(results).transpose()[kk,1])
    y_lower_p[kk] = np.mean(np.array(results).transpose()[kk,1]) - np.std(np.array(results).transpose()[kk,1])
    y_p_mean[kk] = np.mean(np.array(results).transpose()[kk,1])
    #y_upper_qp[kk] = np.mean(var_test_qp[kk, :]) + np.std(var_test_qp[kk, :])
    #y_lower_qp[kk] = np.mean(var_test_qp[kk, :]) - np.std(var_test_qp[kk, :])
    #y_qp_mean[kk] = np.mean(var_test_qp[kk, :]

# Plot using verification process and raw data vs. the theoretical value
plt.figure()
plt.fill_between(x, y_upper, y_lower, color=[0.5, 0.1, 1], alpha=0.3)
plt.plot(x, (2 * V11_aux) * np.ones(len(resol)), linewidth=2, linestyle='--', color=[1, 0.411764705882353, 0.16078431372549])
plt.fill_between(x, y_upper_p, y_lower_p, color=[0.5, 0.1, 1], alpha=0.3)
plt.plot(x, (2 * V22_aux) * np.ones(len(resol)), linewidth=2, linestyle='--', color=[1, 0.411764705882353, 0.16078431372549])
plt.grid(True)
plt.legend(['Experiment', 'Theory'], loc='upper left')
#plt.xticks(ticks=np.arange(1, 121, 10), labels=[str(i) for i in range(1, 121, 10)])
plt.xticks(ticks=[1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120], labels=['1','10','20','30','40','50','60','70','80','90','100','110','120'])
plt.xlim([1, 120])
plt.ylim([0, 3.2e4])
plt.ylabel('Verification Variance')
plt.xlabel('Frequency resolution (Hz)')
plt.show()
