import numpy as np
from math import *
from cmath import *
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.ticker import EngFormatter

# set latex font as default for plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)

def numpify(var, dim = 1):
    """ 
    numpify transforms a list or an integer in a numpy array and returns it

    INPUT
        - var: the list or an integer to be tranformed in a numpy array
        - dim (opt): if var is an integer, the user can specify the dimension of the output array, whose values will be all set to var
    
    OUTPUT
        - if var if a list, the function returns it transformed in a numpy array
        - if var is an integer a no dim is specified, the function returns a single-dimensional array with var
        - if var is an integer a dim is specified, then the function returns a dim-dimensional array whose values are all set to var
    """

    if(isinstance(var, list)):
        var = np.array(var)
    elif(isinstance(var, float) or isinstance(var, int)):
        var = var*np.ones(dim)
    return var

def linreg(x, y, dy=[], dx=[], logx=False, logy=False):
    """ 
    LINREG returns the parameters for a linear fit of the form  y = m*x + b, based on least squares formulas.

    INPUTS:
        x list or array of independent physical quantity
        y list or array of dependent physical quantity
        (opt) dy scalar or array of uncertainties of y
        (opt) dx scalar or array of uncertainties of x
        (opt) logy boolean flag to make the fit with log(y) instead of y
        (opt) logx boolean flag to make the fit with log(x) instead of x
        
    OUTPUT: A dictionary with the following fields
        "m" estimated slope of the line y = m*x + b
        "b" estimated intercept of the line y = m*x + b
        "dm" uncertainty of m
        "db" uncertainty of b
        "chi2r" reduced chi square, i.e. the overall discrepancy between the best straight line and the experimental points
        "dof"  number of degrees of freedom
    
    NOTES
        - x and y must have equal size
        - If no y errors are given, equal errors aìfor all points are assumed and parameters are calculated based on the      assumption chi2 = 1
        - If dx is given, the fit parameters are first computed ignoring dx errors, then the slope is used to convert the
          dx error in an dy one, dy is updated accordingly and the fit is
          performed again
        - If dy or dx are scalars, it is assumed they are equal for all points
    """
    # if x, y, dy, dx are lists, trasform them in numpy arrays
    var = [x, y, dx, dy]
    for i, val in enumerate(var):
        if (isinstance(val, list)):
            var[i] = np.array(var[i])
        elif (isinstance(val, int)):
            var[i] = np.array([var[i]])
    [x, y, dx, dy] = var

    # check that x and y vectors are of equal shape
    if (len(x)!=len(y)):
        print("Error in linreg: x and y have different shapes")
        return
        
    # make vectors if dy or dx are single values
    if (len(dy) == 1):
        dy = dy*np.ones(len(y))
    if (len(dx == 1)):
        dx = dx*np.ones(len(y))

    # make logarithms if specified (and modify accordingly the uncertainties)
    if (logy):
        try:
            y = log(y)
            dy = dy/y
        except:
            print("Could not perform log(y). Check y values")
    if (logx):
        try:
            x = log(x)
            dx = dx/x
        except:
            print("Could not perform log(x). Check x values")
        
    # compute number of points
    N = len(x)

    # attribute weights 
    noerrors = False #default
    if (len(dy)==0):
        noerrors = True
        wt = np.ones(len(y))
    else:
        wt = 1/(dy**2)
    
    # Compose the ingredients for the textbook LSQ solution for linear model
    Sw  = sum(wt)
    Sx  = sum(x * wt)
    Sx2 = sum(x**2 * wt)
    Sy  = sum(y * wt)
    sxy = sum(x * y * wt)
    
    Delta = Sw*Sx2 - Sx**2
  
    # Calculate the linear fit coefficients using the textbook LSQ solution for linear model
    m = (Sw * sxy - Sx * Sy) / Delta
    b = (Sx2 * Sy - Sx * sxy) / Delta
    
    # Calculate the number of degrees of freedom
    # -2 is for 2 degrees of freedom removed by fit, m and b
    dof = N - 2

    if noerrors:
        # The user did not provide dy errors
        # assign the reduced chi^2 a value of 1
        chi2 = 1
        # estimate 'a posteriori' the data uncertainty
        dy_post = sqrt(sum(1 / dof * (y - (m*x + b))**2)) * np.ones(len(y))
        # call the routine using these estimated uncertainties
        out = linreg(x, y, dy=dy_post)
        m, b, dm, db, chi2r, dof = out["m"], out["b"], out["dm"], out["db"], out["chi2r"], out["dof"]    
    else:
        # The user did provide dy errors
        # check correct shape of dy
        if(len(y)!=len(dy)):
            print("Error: dy and y have different shapes")
            return
        dm = sqrt(Sw / Delta)
        db = sqrt(Sx2 / Delta)
        
        # estimate the reduced chi^2
        chi2r = sum((y - (m*x + b))**2 * wt) / dof
        
        # estimate the contribution of the uncertainties along x, if any
        if (len(dx) != 0):
            # check correct shape of dx
            if(len(dx)!= len(y)):
                print("Error: dx and the other variables have different shapes")
                return
            # propagate the dx errors via the first parameter estimates
            dy_est = sqrt(dy^2 + (m * dx)**2)
            
            # call the routine using these updated uncertainties
            out = linreg(x, y, dy=dy_est)
            m, b, dm, db, chi2r, dof = out["m"], out["b"], out["dm"], out["db"], out["chi2r"], out["dof"]

    return {"m":m, "b":b, "dm":dm, "db":db, "chi2r":chi2r, "dof":dof}

def bodeplot(f, H=[], amp=[], Phase=[], figure=[], deg=False, asline=False, plotDeg = True):
    """ 
    BODEPLOT plots the amplitude and phase diagrams of the transfer function given as input
    
    INPUT: The tranfer function can be passed as input in two different ways
        - either as a vector of complex numbers passed as "H = vector"
        - or as two separate vectors of real number indicating the amplitude and the phase, passed as "amp = ampvector, Phase = phasevector"
        Optional input:
        - figure: a figure object to which the lines should be added
        - deg: if True, it is assumed that the given phase is in degrees
        - asline: set it True you want the function to appear as a smooth line. Default is False and isolated points are displayed
        - plotDeg: plots the phase in degrees. Default is True
    
    OUTPUT
        The function creates and returns a matplotlib figure containing two subplots, the first one for the amplitude and the second one for the phase. A logarithmic (base 10) scale is used on the x axis.
        IMPORTANT: The user must use the "show" function of pyplot in matplotlib to display the figure.
    """

    # if phase id given in degrees, transform it in radians
    if(deg):
        Phase = Phase * pi/180

    # calculate modulus and phase of the transfer function if complex H is given
    if (len(H) != 0):
        amp = abs(H)
        Phase = np.angle(H)

    #check correct shapes
    if(len(f) != len(amp)):
        print("Error: frequence and amp have different shapes")
        return
    if(len(f)!=len(Phase)):
        print("Error: frequence and phase have different shapes")
        return

    # if the figure parameter is given, add the plot to that figure, otherwise create a figure
    if(figure==[]):
            figure,_ = plt.subplots(nrows=2, ncols=1)
    [ampax, phaseax] = figure.get_axes()

    # amplitude plot
    ampax.set_xscale("log")
    ampplot = ampax.plot(f, amp)
    ampax.grid(b=True, which="both")       # no idea what "b=True" does, but without it the grid doesn't show up

    # amplitude style setup    
    if(not asline):
        plt.setp(ampplot, ls = ' ', c = "red", marker='o', ms=4)
    else:
        plt.setp(ampplot, ls = '-', c = "black")
    ampax.set_xlabel(r"$f$ [Hz]")
    ampax.set_ylabel(r"$|H|$")


    # phase plot
    phaseax.set_xscale("log")
    if(plotDeg):
        Phase = Phase*180/pi
    phaseplot = phaseax.plot(f, Phase, color="black")
    phaseax.grid(b=True, which="both") 

    # phase style setup
    if(not asline):
        plt.setp(phaseplot, ls = ' ', c = "red", marker='o', ms=4)
    else:
        plt.setp(phaseplot, ls = '-', c = "black")
    phaseax.set_xlabel(r"$f$ [Hz]")
    phaseax.set_ylabel(r"$\phi$")
    if(plotDeg):
        phaseax.yaxis.set_major_formatter(EngFormatter(unit=u"°"))

    # set distance between subfigures
    plt.subplots_adjust(hspace = .5)
    return figure


