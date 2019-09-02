import numpy as np
import sys
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------#

'''
This program should, given the proper parameters provide the information to
select a filter and exposure time in addition to target selection.

To Do:
- change to a kwarg system instead of raw input
- and help function 
'''

#------------------------------------------------------------------------------#

# Plank function 

def Blam( Lambda, T):
    h = np.float64(6.62607004e-34)
    c = np.float64(299792458)
    k = np.float64(1.38064852e-23)
    bot = np.expm1((h*c)/(np.float64(Lambda)*k*np.float64(T)))
    top = 2*h*c**2*np.pi
    first = top/(np.float64(Lambda)**5)
    second = 1/bot
    Bl = first*second 
    return(Bl)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Generating Plank Curves

# wavelengths 
lam = np.linspace(5.89e-7,7.27e-7,10000) # red filter 
lamf = np.linspace(0.5e-9,1.0e-6,10000) # full spectrum
lamv = np.linspace(5.07e-7,5.95e-7,10000) # green filter
lamb = np.linspace(3.98e-7,4.92e-7,10000) # blue filter

T = float(input('Temperature of star in Kelvins: ') )
m = float(input('Apparent AB magnitude of Star; ') ) 

# Running Plank Function 
Bl = Blam(lam,T)
Blf = Blam(lamf,T)
Blv = Blam(lamv,T)
Blb = Blam(lamb,T)

# Ploting the Plank Curve
plt.figure()
plt.plot(lamf*10**(9),Blf*10**(-9),'--',color='black')
plt.plot(lamv*10**(9),Blv*10**(-9),color='green',label='Green Filter')
plt.fill_between(lamv*10**(9),Blv*10**(-9),0,color='green',alpha=0.4)
plt.plot(lam*10**(9),Bl*10**(-9),color='red',label='Red Filter')
plt.fill_between(lam*10**(9),Bl*10**(-9),0,color='red',alpha=0.4)
plt.plot(lamb*10**(9),Blb*10**(-9),color='blue',label='Blue Filter')
plt.fill_between(lamb*10**(9),Blb*10**(-9),0,color='blue',alpha=0.4)
plt.ylim(0)
plt.xlim(0,1000)
plt.legend()
plt.xlabel('Wavelength (nm)',fontsize= 10) 
plt.ylabel('Spectral Radiance ($ W m^{-2}  nm^{-1}$)', fontsize = 10)
plt.title('%.0fk Theoretical Spectral Black Body Profile' %( T), fontsize = 15)
plt.show()

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Numerical Integration

h = np.float64(6.62607004e-34)
c = np.float64(299792458)

Lambda = [lam,lamb,lamv,lamf]
Blrange = [Bl,Blb,Blv,Blf]
F = [0 , 0, 0, 0]

for i in range(len(Lambda)):

    a = []
    stepsize = max(Lambda[i])-min(Lambda[i])
    
    for j in range(len(Lambda[i])):
        z = Blrange[i]
        area = stepsize*z[j]
        a.append(area)
        
    F[i] = sum(a)

Fr = F[0]
Fb = F[1]
Fv = F[2]
Ffull = F[3]

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Filter Scaling

x = input('Please input desired filter r,v,b:' )

if x == 'r' :
    Lratio = Fr/Ffull
    Ep = h*c/(6.58e-6)
elif x == 'v' :
    Lratio = Fv/Ffull
    Ep = h*c/(5.51e-6)
elif x == 'b' :
    Lratio = Fb/Ffull
    Ep = h*c/(4.45e-6)

appF = 10**((m+48.6)/(-2.5))

FilterF = appF*Lratio
nps = FilterF/Ep
print(nps)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

 


 









