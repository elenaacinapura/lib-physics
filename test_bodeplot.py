from libphysics import *

freq = numpify([30, 50, 80, 150, 400, 700, 1e3, 2e3, 3e3,   5e3, 8e3, 12e3, 17e3, 20e3])
Vin = numpify([1.88, 1.87, 1.87, 1.87, 1.87, 1.87, 1.86, 1.85, 1.84,   1.84, 1.83, 1.83, 1.83, 1.83])
Vout = numpify([1.84, 1.84, 1.84, 1.83, 1.78, 1.69, 1.62, 1.1, .838, .556, .374, .265,.206, .182])
fase = np.array([1.4, 2, 3, 5, 14, 23.5, 31, 48.7, 58,  63, 66, 64, 60, 58 - 360])*-1

modH = Vout/Vin
b1 = bodeplot(freq, Amp=modH, Phase=fase, deg=True, color = "blue", asline=True, linestyle='-.', linear_yscale=1)
plt.show()