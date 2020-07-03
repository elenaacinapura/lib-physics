from libphysics import *

freq = [30, 50, 80, 150, 400, 700, 1e3, 2e3, 3e3,   5e3, 8e3, 12e3, 17e3, 20e3]
Vin = [1.88, 1.87, 1.87, 1.87, 1.87, 1.87, 1.86, 1.85, 1.84,   1.84, 1.83, 1.83, 1.83, 1.83]
Vout = [1.84, 1.84, 1.84, 1.83, 1.78, 1.69, 1.62, 1.1, .838, .556, .374, .265,.206, .182]
fase = np.array([1.4, 2, 3, 5, 14, 23.5, 31, 48.7, 58,  63, 66, 64, 60, 58])*-1


Vin = numpify(Vin)
#Vout = numpify(Vout)
fase = numpify(fase)
freq = numpify(freq)
modH = Vout/Vin


f = np.hstack([numpify(Vin*3, column=True), numpify(freq*4, column= True)])

out = lsq_fit(Vout, f, fase)
print(out["fit_out"])