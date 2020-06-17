import numpy as np
import csv
import os
from math import *
from cmath import *
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.ticker import EngFormatter

# set latex font as default for plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)

def numpify(var, dim = 1, column=True):
    """ 
    numpify transforms a list or an integer in a numpy ndarray column and returns it

    INPUT
        - var: the list or an integer to be tranformed in a numpy array
        - dim (opt): if var is an integer, the user can specify the dimension of the output array, whose values will be all set to var
    
    OUTPUT
        - if var if a list, the function returns it transformed in a numpy array
        - if var is an integer a no dim is specified, the function returns a single-dimensional array with var
        - if var is an integer a dim is specified, then the function returns a dim-dimensional array whose values are all set to var
    """

    if(isinstance(var, list)):
        if(column):
            var = np.array([var]).transpose()
        else: 
            var = np.array(var)
    elif(isinstance(var, float) or isinstance(var, int)):
        if(column):
            var = var*np.ones((dim,1))
        else:
            var = var*np.ones(dim)
    elif(column and len(np.shape(var))==1):
            var = var.reshape(len(var),1)
    return var

def readCSV(file, skiprows=0, cols=[], untilrow=0):
    """
    readCSV reads the content of a csv or text file and returns its columns as arrays

    INPUT
        - file: string containing the name of the file to be read (relative path)
        - skiprows (opt): number of lines to be skipped at the beginning of the file (e.g. because there's a header with names for the columns). By default every row is read
        - cols (opt): array with the indexes of the columns to be read (start counting with index zero please). By default every column is read

    OUTPUT
        - a list of arrays, each one containing the content of a column of the file
    """
    # open the input file
    filetoread = os.path.join(file)
    if os.path.isfile(filetoread):
        with open(file, 'r') as f:
            reader = csv.reader(f)

            # count number of columns if not given ho many to count
            if (cols==[]):
                ncols = len(next(reader)) # Read first line and count columns
                cols = [i for i in range(ncols)]
            else:
                ncols = len(cols)   
            # return to the beginning of the file
            f.seek(0) 

            # data structure to store the input
            data = np.ndarray((1, ncols))

            # loop on the lines of the file skipping rows if told so
            for i,row in enumerate(reader):
                if (i<skiprows):
                    continue
                if (untilrow != 0 and i>= untilrow):
                    break
                # make a list from the line (reading only the wanted columns)
                r = []
                for j, element in enumerate(row):
                    if(j in cols):
                        try:
                            r.append(float(element))
                        except:
                            print("Couldn't read input in row ", i, ", column ", j)
                            continue
                if (i==0+skiprows):
                    data[0] = r
                else:
                    try:
                        data = np.vstack([data, r])  
                    except:
                        continue                    
    else:
        print("Error: couldn't find file " + file + ". Make sure to execute this script in the same folder of the file to read")
        return
    
    # return a list of separate columns
    output = []
    for i in range(ncols):
        output.append(data[:,i])
    
    return output

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

def bodeplot(f, H=[], Amp=[], Phase=[], figure=[], deg=True, err=False, Amperr=[], Phaseerr=[],asline=False, plotDeg = True):
    """ 
    BODEPLOT plots the amplitude and phase diagrams of the transfer function given as input
    
    INPUT: The tranfer function can be passed as input in two different ways
        - either as a vector of complex numbers passed as "H = vector"
        - or as two separate vectors of real number indicating the amplitude and the phase, passed as "Amp = ampvector, Phase = phasevector"
        Optional input:
        - figure: a figure object to which the lines should be added
        - deg: if True, it is assumed that the given phase is in degrees
        - asline: set it True you want the function to appear as a smooth line. Default is False and isolated points are displayed
        - plotDeg: plots the phase in degrees. Default is True
    
    OUTPUT
        The function creates and returns a matplotlib figure containing two subplots, the first one for the amplitude and the second one for the phase. A logarithmic (base 10) scale is used on the x axis.
        IMPORTANT: The user must use the "show" function of pyplot in matplotlib to display the figure.
        To work with the axes of the figure, type [ax1,ax2] = figure.axes
    """

    # if phase is given in degrees, transform it in radians
    if(Phase!=[] and deg):
        Phase = Phase * pi/180

    # calculate modulus and phase of the transfer function if complex H is given
    if (len(H) != 0):
        Amp = abs(H)
        Phase = np.angle(H)

    #check correct shapes
    if(len(f) != len(Amp)):
        print("Error: frequence and Amp have different shapes")
        return
    if(len(f)!=len(Phase)):
        print("Error: frequence and phase have different shapes")
        return
    if(err):
        if(len(f)!=len(Amperr)):
            print("Error: frequence and amp error have different shapes")
            return
        if(len(f)!=len(Phaseerr)):
            print("Error: frequence and phase error have different shapes")
            return

    # if the figure parameter is given, add the plot to that figure, otherwise create a figure
    if(figure==[]):
            figure,_ = plt.subplots(nrows=2, ncols=1)
    [ampax, phaseax] = figure.get_axes()
    

    # amplitude plot
    ampax.set_xscale("log")
    if(err):
        ampplot = ampax.errorbar(f, Amp, yerr=Amperr, ls=' ', c="red", marker="o", ms=4, ecolor="red")
    else: 
        ampplot = ampax.plot(f, Amp)
    ampax.grid(b=True, which="both")       # no idea what "b=True" does, but without it the grid doesn't show up

    # amplitude style setup    
    if(not asline and not err):
        plt.setp(ampplot, ls = ' ', c = "red", marker='o', ms=4)
    elif(asline):
        plt.setp(ampplot, ls = '-', c = "black")
    ampax.set_xlabel(r"$f$ [Hz]")
    ampax.set_ylabel(r"$|H|$")


    # phase plot
    phaseax.set_xscale("log")
    if(plotDeg):
        Phase = Phase*180/pi
    if(err):
        phaseplot = phaseax.errorbar(f, Phase, yerr=Phaseerr, ls = ' ', c = "red", marker='o', ms=4)
    else: 
        phaseplot = phaseax.plot(f, Phase, color="black")
    phaseax.grid(b=True, which="both") 

    # phase style setup
    if(not asline and not err):
        plt.setp(phaseplot, ls = ' ', c = "red", marker='o', ms=4)
    elif(asline):
        plt.setp(phaseplot, ls = '-', c = "black")
    phaseax.set_xlabel(r"$f$ [Hz]")
    phaseax.set_ylabel(r"$\phi$")
    if(plotDeg):
        phaseax.yaxis.set_major_formatter(EngFormatter(unit=u"°"))

    # set distance between subfigures
    plt.subplots_adjust(hspace = .5)
    return figure

def lsq_fit(y, f, dy):
    """ 
    lsq_fit fits data y as a sum of functions f_j faccording to the model y = SUM_j(a_j * f_j)
    
    INPUT:
        y: array of dependent data
        f: matrix where the j-th column is the expected value for f_j(x), the i-th row represent the i-th component
        dy: array of uncertainties for y
    
    OUTPUT: the vector of the coefficients a_j that best fit y = SUM_j(a_j * f_j)
    """
    # check correct shapes
    if (len(y) != len(f)):
        print("Error: the size of y doesn't equal the number of rows of f")
        return
    y = numpify(y)
    dy = numpify(dy)
    # print(y)
    # print(dy)
    n_points = len(y)
    n_func = np.shape(f)[1]

    V = np.ndarray((n_func,1))
    G = np.ndarray((n_func, n_func))

    for i in range(n_func):
        #print(f[:,[i]]* y / (dy)**2)
        V[i] = np.sum(f[:,[i]]* y / (dy)**2)
        for j in range(n_func):
            G[i,j] = np.sum(f[:,[i]] * f[:,[j]] / (dy)**2)
    #print(V)
    C = np.linalg.inv(G)
    fit_out = C.dot(V)
    dfit_out = np.ndarray(np.shape(fit_out))
    for i in range(n_func):
        dfit_out[i] = C[i,i]


    # calculate discrepancy
    y_fit = f.dot(fit_out)
    if (np.shape(y) != np.shape(y_fit)):
        y = y.reshape(np.shape(y_fit))
    if(np.shape(dy) != np.shape(y_fit)):
        dy = dy.reshape(np.shape(y_fit))
    dy_res = y - y_fit
    n_dof = len(y) - len(fit_out)
    chi2r = np.sum(dy_res**2 / dy**2) / n_dof

    return {"fit_out":fit_out, "dfit_out":dfit_out, "chi2r":chi2r}

