# via_metabolica
from cinetica_secondo_ordine import *
from cinetica_serie_2_1 import *

t_inizio=0
t_fine=4
t_campioni=5000
time=np.linspace(t_inizio,t_fine,t_campioni)
NO=5.84e-7
k1=8e7
sCG=2e-6
k2=1.3

#Esempio di reazione chimica con cinetica del secondo ordine, A+B->Prodotto
#Grafici delle concentrazione dei substrati e del prodotto
decadimento_NO=concentrazione_A(NO,sCG,k1,time)
decadimento_sCG=concentraziione_B(NO,sCG,k1,time)
prodotto_cinetica_secondo_ordine=concentrazione_P(NO,sCG,k1,time)

plt.subplot(3,1,1)
plt.plot(time,decadimento_NO)
plt.subplot(3,1,2)
plt.plot(time,decadimento_sCG)
plt.subplot(3,1,3)
plt.plot(time,prodotto_cinetica_secondo_ordine)
plt.show()

#Invece noi stiamo considerando due reazioni chimiche in serie. La prima reazione ha una cinetica del secondo ordine, mentre la seconda ha una cinetica del primo ordine
#A+B->I->P
#Questa funzione vuole rappresentrare la concentrazione di I nel tempo
#Si ricava risolvendo l'equazione differenziale dI/dt=-k2*I+k1*A*B quindi una equazione differenziale lineare del primo ordine

concentrazione_prodotto_intermedio_I=prodotto_intermedio_serie_2_1(NO,sCG,k1,k2,t_inizio,t_fine,0,'RK45',time)
plt.plot(time, concentrazione_prodotto_intermedio_I)
plt.show()

#Per il bilancio di massa vale: P=A0+B0-A-B-I

concentrazione_prodotto_P=np.zeros(t_campioni)
for i in range(t_campioni):
    concentrazione_prodotto_P[i]=NO+sCG-decadimento_NO[i]-decadimento_sCG[i]-concentrazione_prodotto_intermedio_I[i]

plt.plot(time,concentrazione_prodotto_P)
plt.show()
