#Import Libarys
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
from scipy import optimize
from numpy import sqrt, sin, cos, pi

#Programm to analyze data from Sinusoidal testing, based on the work of Todd C. Rasmussen, Kevin G. Haborak and Michael H. Young
#Written by Golo Voelker at the University of Tübingen to simplyfy and automize the analytical estimation of Transmissivity, Storativity and Diffusivity from Sinusoidal
#pumping test data
#Needes Libarys: Numpy, Matplotlib,scipy
#Version 1.0
#Date 02.08.2021
#No guarantee for completeness or correctness
#Help in the case that phi0 is negative: Delete some of the inputdata, might sound rather strange, but results in a change in phi0


#Input of known constances (Periodlength and r)
print("The inputfiles containing Q and s have to be named Q.txt and S.txt. They needt to be within the same folder")
p = int(input("Enter Periodlength in seconds: "))
r = int(input("Distance from Source to observation well in m: "))             #Distance of the wells
w = (2*pi)/p                                                        #frequency of the pumping


#The input is built around the output from the Software nSIGHTS, meaning time in coloumn 0 is in excel days, the Pumpingrate can have any unit
#The input is deleting the first line to be optimized using the output from nSIGHTS
#Make sure to check the fit for irregularities
#Get Q1 and Q2

rawQ = np.loadtxt('Q.txt',skiprows=1)
tq = rawQ[:,0]
for i in range (1,tq.shape[0]):      #set start time to 0
    tq[i] = tq[i]-tq[0]
tq[0] = 0
tq = tq*86400                         #convert from exceltime to seconds
Q = rawQ[:,1]
   
def func(tq, q1, q2):                #calculation of q1 and q2 using least linear squres as described in equation 26
    return q1*np.cos(w*tq)+q2*np.sin(w*tq)


q1_guess = np.random.rand(1)
q2_guess = np.random.rand(1)    
params, cov = optimize.curve_fit(func, tq, Q, [q1_guess,q2_guess])#Fitting

q1 = params[0]
q2 = params[1]
Q_calc = q1*np.cos(w*tq)+q2*np.sin(w*tq)#calculate fitted function



#The input is built around the output from the Software nSIGHTS, meaning time in coloumn 0 is in excel days, the Drawdown is in kPa
#The input pressure is corrected for the offset and converted from kPa to m watercoloumn
#The input is deleting the first line to be optimized using the output from nSIGHTS
#Make sure to check the fit for irregularities
#Get s1 and s2

rawS = np.loadtxt('S.txt',skiprows=1)
ts = rawS[:,0]
for i in range (1,ts.shape[0]):      #set start time to 0
    ts[i] = ts[i]-ts[0]
ts[0] = 0
ts = ts*86400                         #convert from exeltime to seconds
S = rawS[:,1]
#transfer pressure to drawdown [m] and correct offset
meanS = np.mean(S)
S = (S - meanS)/10                  #convert from kPa to m of watercoloumn

def func(ts, s1, s2):                #calculation of s1 and s2 using least linear squres as described in equation 29
    return s1*np.cos(w*ts)+s2*np.sin(w*ts)
s1_guess = np.random.rand(1)
s2_guess = np.random.rand(1)    
params, cov = optimize.curve_fit(func, ts, S, [s1_guess,s2_guess])#Fitting

s1 = params[0]
s2 = params[1]
S_calc = s1*np.cos(w*ts)+s2*np.sin(w*ts)#calculate fitted function


#Plot Drawdown and Pumpingrate with corresponding fits

print("close figure to proceed")

plt.figure(1)
plt.subplot(211)
plt.plot(ts,S, label="drawdown")
plt.plot(ts,S_calc, label="s_drawdown")
plt.xlabel('time [s]')
plt.ylabel('Drawdown [m]')
plt.title('Fitting of the Drawdown')

plt.subplot(212)
plt.plot(tq,Q, label="Pumpingrate")
plt.plot(tq,Q_calc, label="s_Pumpingrate")
plt.xlabel('time [s]')
plt.ylabel('Pumpingrate [m^3/s]')
plt.title('Fitting of the Pumpingrate')
plt.show()


#Calculate s0, Q0 and phi0
Q0 = np.sqrt(q1**2+q2**2)
s0 = np.sqrt(s1**2+s2**2)/Q0
phi0 = np.arctan(s2/s1)-np.arctan(q2/q1)

#calculate u
c = np.array([-0.12665,2.8642,-0.47779,0.16586,-0.076402,0.03089])
lnu_old = 0
for i in range (0,6):
    ci = c[i]
    lnu_new = c[i]*np.log(phi0)**i
    lnu_old = lnu_old+lnu_new
print(lnu_old)
u = np.exp(lnu_old)

#Calculate Diffusivity
D=(w*r**2)/u
print("Diffusivity [m^2/s]:")
print(D)
#calculate the Transmissivity
T = special.k0(np.sqrt(np.absolute(1j+u)))/(2*pi*s0)
print("Transmissivity [m^2/s]:")
print(T)
#calculate S
S = T/D
print("Storativity [-]:")
print(S)
print("")
print("")
input("Press enter to exit")