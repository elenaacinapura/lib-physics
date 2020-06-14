from libphysics import *

# freq = [30, 50, 80, 150, 400, 700, 1e3, 2e3, 3e3,   5e3, 8e3, 12e3, 17e3, 20e3]
# Vin = [1.88, 1.87, 1.87, 1.87, 1.87, 1.87, 1.86, 1.85, 1.84,   1.84, 1.83, 1.83, 1.83, 1.83]
# Vout = [1.84, 1.84, 1.84, 1.83, 1.78, 1.69, 1.62, 1.1, .838, .556, .374, .265,.206, .182]
# fase = np.array([1.4, 2, 3, 5, 14, 23.5, 31, 48.7, 58,  63, 66, 64, 60, 58])*-1


# Vin = numpify(Vin)
# Vout = numpify(Vout)
# fase = numpify(fase)
# freq = numpify(freq)
# modH = Vout/Vin

# f = np.hstack([Vin*3, freq*4])

# out = lsq_fit(Vout, f, fase)
# print(out)

# b1 = bodeplot(freq, amp=modH, Phase=fase, deg=True)

# C = 33.6e-9
# R1 = (.98577+ 1.0022 + .9969)*1e3
# R2 = .11933e3

# H_teo = (1+1j*2*pi*freq*C*R2)/(1 + 1j*2*pi*freq*C*(R1 + R2))

# b2 = bodeplot(freq, H=H_teo, asline=True, figure=b1)

# [x,y,z] = readCSV("input_test.txt")
# print(x, y, z)


[f, amp, fase, damp, dfase] = readCSV("input_test.csv", skiprows=1)
bodeplot(f, Amp=amp, Phase=fase, err=True, Amperr=damp, Phaseerr=dfase)
plt.show()