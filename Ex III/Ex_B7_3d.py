from scipy import stats
import numpy as np
from matplotlib import pyplot as plt

x = np.linspace(-0.25,0.25,500)#in meter
#x = np.sqrt(1.602e-19)*x
m = 1.67372e-27# kg Masse H-Atom
h_qu = 4.135667e-15/2/np.pi #eV*s H-quer
omega = 8.5e13# Hz kreisfrequenz

alph = np.sqrt(m*omega/h_qu)
anf= np.sqrt(alph)/(np.pi)**(1/4)
e_const = np.exp(-alph**2*x**2/2)#letzter Teil der Gleichung



H_0 = 1
H_1 = 2*alph*x
H_2 = 4*alph**2*x**2
H_3 = 16*alph**4*x**4-48*alph**2*x**2+12

E_0 = h_qu*omega*0.5
E_1 = h_qu*omega*3/2
E_2 = h_qu*omega*5/2
E_3 = h_qu*omega*7/2

E_0_plot = np.zeros(shape=(len(x)))
for i in range(len(E_0_plot)):
    E_0_plot[i] = E_0

E_1_plot = np.zeros(shape=(len(x)))
#print(E_1_plot)
for i in range(len(E_1_plot)):
    E_1_plot[i] = E_1
#print(E_1_plot)

E_2_plot = np.zeros(shape=(len(x)))
#print(E_2_plot)
for i in range(len(E_2_plot)):
    E_2_plot[i] = E_2
#print(E_2_plot)

E_3_plot = np.zeros(shape=(len(x)))
#print(E_3_plot)
for i in range(len(E_3_plot)):
    E_3_plot[i] = E_3
#print(E_3_plot)

psi_0 = anf*1*H_0*e_const*0.005 + E_0_plot#factor 0.005 to get well relationship between both plots
Psi_1 = anf*1/(np.sqrt(2))*H_1*e_const*0.005 + E_1_plot
Psi_2 = anf*1/(np.sqrt(8))*H_2*e_const*0.005 + E_2_plot
Psi_3 = anf*1/np.sqrt(2**3*6)*H_3*e_const*0.004 + E_3_plot

potential = 0.5*x**2*omega**2*m#nach Aufgabe c) 

for i in range(len(x)):
    print(Psi_1[i])

#print(x)

plt.plot(x, psi_0)
plt.plot(x, Psi_1)
plt.plot(x, Psi_2)
plt.plot(x, Psi_3)
plt.plot(x,potential)
plt.plot(x, E_0_plot, label = "E=0")
plt.plot(x, E_1_plot, label = "E=1")
plt.plot(x, E_2_plot, label = "E=2")
plt.plot(x, E_3_plot, label = "E=3")
plt.xlabel("x in nm")
plt.xlim(-0.25,0.25)
plt.ylabel("eV")


#plt.plot(x, y, label = "Minima", color= 'green', linewidth=1)


#plt.title("7.3d))

plt.savefig("7_3d_Molekuelschwingung.pdf", format="PDF")
#plt.show()
##plt.close()

####Number 7_1
q = 1.602e-19#C
m = 0.3* 9.109e-31#kg
#phi= 4*1.602e-19#J
phi = 4#V
d=4.1e-9#nm
h_quer= 6.626e-34/2/np.pi#Js
U= np.array([5,10,15])
T = 300#K
kb = 1.38e-23#J/K

P = np.exp(-4/3*np.sqrt(2*q*m)*phi**(3/2)*d/h_quer/U)
print(P)

P_0 = np.exp(-2*d/h_quer*np.sqrt(2*m*q*phi))

print("P_0",P)
Tau_0 = 2*4.1e-9/8.065e17
print("Tau_0",Tau_0)
Tau_0V = Tau_0/P_0
Tau = Tau_0/P
print("Lebenszeit",Tau_0V,Tau)
print("Zeit zu der 0.1%, 0.5% der Elektronen noch da sind = ", -np.log(0.1)*Tau_0V,-np.log(0.5)*Tau_0V)
print("Zeit zu der 0.1% der Elektronen noch da sind = ,wenn Spannung angelegt wird", -np.log(0.1)*Tau)