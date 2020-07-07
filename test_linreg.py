from libphysics import *

x = numpify([1.1,2,3,4,5,6,7,8])
dy = np.ones(len(x))*0.1
dx = np.ones(len(x))* 0.2
y = 3*x + 2
# otherwise fit is too perfect
y[0] += 0.1
y[1] += 0.2

fit = linreg(x, y, dy = dy, dx = dx)
print(fit["m"])
print(fit["b"])