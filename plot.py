
# coding: utf-8

# In[1]:


from math import*
from scipy import*
import numpy as np
import matplotlib.pyplot as plt




# In[2]:


fichier = np.loadtxt("Vec_psi.txt")
fichier2 = np.loadtxt("Vec_N.txt")
fichier3 = np.loadtxt("Vec_P.txt")
fichier4 = np.loadtxt("Vec_psi_eq.txt")
fichier5 = np.loadtxt("Vec_n_eq.txt")
fichier6 = np.loadtxt("Vec_p_eq.txt")


# In[3]:


plt.plot(fichier[:,0],fichier[:,1],)
plt.title("Vecteur psi à t=0.1")
plt.xlabel("x")
plt.ylabel("psi")
plt.show()
plt.plot(fichier2[:,0],fichier2[:,1])
plt.title("Vecteur n à t=0.1")
plt.xlabel("x")
plt.ylabel("n")
plt.show()
plt.plot(fichier3[:,0],fichier3[:,1])
plt.title("Vecteur p à t=0.1")
plt.xlabel("x")
plt.ylabel("p")
plt.show()
plt.plot(fichier[:,0],fichier4[:,1],)
plt.title("Vecteur psi à l'equilibre")
plt.xlabel("x")
plt.ylabel("psi")
plt.show()
plt.plot(fichier2[:,0],fichier5[:,1])
plt.title("Vecteur n à l'equilibre")
plt.xlabel("x")
plt.ylabel("n")
plt.show()
plt.plot(fichier3[:,0],fichier6[:,1])
plt.title("Vecteur p à l'équilibre")
plt.xlabel("x")
plt.ylabel("p")
plt.show()


# In[9]:


t=[0.00001,0.0005 ,0.001,0.008,0.02,0.04, 0.06 , 0.08 ,0.1 ,0.25   ,0.5 ]
E=[0.46 ,0.45   ,0.44 ,0.373,0.28,0.18, 0.11 ,0.078 ,0.05,0.003,0.00003]
I=[43.14 ,16.75 ,14.7,9.7,7.16,4.6, 2.9 ,1.9,1.24,0.06,0.0006]
plt.plot(t,I)
plt.title('Dissipation en fonction du temps')
plt.xlabel('temps')
plt.ylabel('dissipation')
plt.show()
plt.plot(t,E)
plt.title('Energie en fonction du temps')
plt.xlabel('temps')
plt.ylabel('Energie')
plt.show()

