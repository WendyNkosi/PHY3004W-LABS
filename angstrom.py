import numpy as np
from scipy.optimize import curve_fit    
import matplotlib.pyplot as plt
import math
from scipy.integrate import simps
from scipy.signal import square

file1=np.loadtxt('FINAL.txt')
t=[]
tp=[]
tq=[]
for line in file1:
    t.append(float(line[0]))
    tp.append(float(line[2]))
    tq.append(float(line[3]))

time=np.array(t)
print(len(time))
tempp=np.array(tp)
tempq=np.array(tq)
Nsamples=len(time)
L= 0.062          #distance between thermocouples
T= 750            #time for one cycle
dt= T/Nsamples   # with five cycles, T/Nsamples    #our real data has r=sqrt(2/3)
p=8.45*10**3      #density
c= 385            #specific heat capacity

results=[]
A0=[]
an=[]
bn=[]
dn=[]
phi=[]
wave_num=[]
att_num=[]
Dn=[]
const=[]
n=0

while n!=4:
    w=(n*2*np.pi)/T
    int_tempp=0
    int_tempq=0
    intb_tempp=0
    intb_tempq=0
    n+=1
    for i in range(Nsamples):
        int_tempq+=tempq[i]*dt*np.cos(w*time[i])
        int_tempp+=tempp[i]*dt*np.cos(w*time[i])
        intb_tempq+=tempq[i]*dt*np.sin(w*time[i])
        intb_tempp+=tempp[i]*dt*np.sin(w*time[i])
        const+= np.sin(w*time[i])
    an.append([2*int_tempp/T,2*int_tempq/T])
    bn.append([2*intb_tempp/T,2*intb_tempq/T])
  
    dn.append([np.sqrt((2*int_tempp/T)**2+ (2*intb_tempp/T)**2),np.sqrt((2*int_tempq/T)**2+ (2*intb_tempq/T)**2)]) 
    phi.append([np.arctan2((2*int_tempp/T),(2*intb_tempp/T)), np.arctan2((2*int_tempq/T),(2*intb_tempq/T))])
    wave_num.append(((np.arctan((2*int_tempp/T)/(2*intb_tempp/T)))-(np.arctan((2*int_tempq/T)/(2*intb_tempq/T))))/L)
    att_num.append(np.log((np.sqrt((2*int_tempp/T)**2+ (2*intb_tempp/T)**2))/(np.sqrt((2*int_tempq/T)**2+ (2*intb_tempq/T)**2)))/L)
    Dn.append(w/(2*(((np.arctan((2*int_tempp/T)/(2*intb_tempp/T)))-(np.arctan((2*int_tempq/T)/(2*intb_tempq/T))))/L)*(np.log((np.sqrt((2*int_tempp/T)**2+ (2*intb_tempp/T)**2))/(np.sqrt((2*int_tempq/T)**2+ (2*intb_tempq/T)**2)))/L)))
                  
##print("AN:",an)
##print("BN:",bn)
##print("DN:",Dn)
##print("dn:",dn)
##print("att_num:",att_num)
##print("wave_num:",wave_num)
print("phase difference:",phi)    
##calculating the conductivity
###propagating uncertainty    
cond =[]
for j in Dn:
    cond.append(j*p*c)
print("conductivity:",cond)
max_temp_P = max(tempp)
max_temp_Q = max(tempq)
index_max_temp_P = np.where(tempp == max_temp_P)[0][0]
index_max_temp_Q = np.where(tempq == max_temp_Q)[0][0]

u_an_list = []
u_bn_list = []
for n in range(1, 5):    
    Oa_P = abs(((n*(2*np.pi/750))**2)*max_temp_P*np.cos((n*(2*np.pi/750))*index_max_temp_P))
    #Oavalue_Q = abs(((n*(2*np.pi/750))**2)*max_temp_Q*np.cos((n*(2*np.pi/750))*index_max_temp_Q))
    Ob_P = abs(((n*(2*np.pi/750))**2)*max_temp_P*np.sin((n*(2*np.pi/750))*index_max_temp_P))
    M_an = (Oa_P*T**2)/(12*Nsamples**2)
    M_bn = (Ob_P*T**2)/(12*Nsamples**2)
    u_an_list.append(M_an)
    u_bn_list.append(M_bn)

    
u_dn_list = [] 
u_alphan_list = []  
for y in range(4):
    unc_dn = dn[y]*np.sqrt((u_an_list[y]/an[y])**2 + (u_bn_list[y]/bn[y])**2)
    unc_alphan = phi[y]*np.sqrt((u_an_list[y]/an[y])**2 + (u_bn_list[y]/bn[y])**2)
    u_alphan_list.append(unc_alphan)
    u_dn_list.append(unc_dn)
    
print('Uncertainties in the dns ', u_dn_list)
print('Uncertainties in the conductivty ',u_alphan_list)

u_Dn_list = []
u_phin_list = []
for kk in range(4):
    unc_Dn = Dn[kk]*np.sqrt((u_dn_list[kk]/dn[kk])**2 + (u_alphan_list[kk]/phi[kk]))
    unc_phin = cond[kk]*np.sqrt((unc_Dn/Dn[kk])**2)
    u_Dn_list.append(unc_Dn)
    u_phin_list.append(unc_phin)
print('Uncertainties in the Dn ',u_Dn_list)
print('Uncertainties in the phi ',u_phin_list)
    















##superimposing dns
ya=[]
yb=[]
for i in range(len(time)):
    a = 0
    b = 0
    T_a = 0
    T_b = 0
    count = 1
    for n in range(len(Dn)):
        w=(n*2*np.pi)/T
        a += (dn[n][0])*np.sin(w*time[i]+phi[n][0])
        b += (dn[n][1])*np.sin(w*time[i]+phi[n][1])
    T_a = an[0][0] + a
    T_b = an[0][1] + b
    ya.append(T_a)
    yb.append(T_b) 
    
#plots
##plt.plot(time, tempp,'.',color='red',label='Q')
##plt.plot(time, tempq,'.',color='blue',label='P')
    
#plt.plot(time, ya, color='black',label='Q')
#plt.plot(time, amplitude, color='black',label='P selected')
##plt.xlabel("Time (in s)")
##plt.ylabel("Temperature (in Celcius)")
##plt.title('Fourier series')
##plt.legend()
##plt.show()


