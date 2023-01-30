#===========================================
#>>>>>>>>>>>>>>>>>BIBLIOTECAS<<<<<<<<<<<<<<<
#===========================================

import numpy as np                     # biblioteca de matrizes
import matplotlib.pyplot as plt        # biblioteca de gráfcos
from scipy.sparse.linalg import eigsh  # biblioteca para Diagonalizar a matriz hamiltoniana
from numba import njit                 # biblioteca de performace
from scipy import constants as cte     # biblioteca de constantes físicas
from scipy.sparse import diags         # biblioteca de cria uma matriz tridiagonal (sparse)

#===========================================
#>>>>>>>>>>>>>>>>>CONSTANTES<<<<<<<<<<<<<<<<
#===========================================

Ry = cte.value(u'Rydberg constant times hc in eV')          # Constante de Rydberg ---> eV
A0 = cte.value(u'Bohr radius')/cte.value(u'Angstrom star')  # Raio de Bohr ---> Angstroms 
k = 3                                                       # Número de autoenergias a serem encontradas
typem = np.float64

#===========================================
#>>>>>>>>>>>>>>>>>>>>GRID<<<<<<<<<<<<<<<<<<<
#===========================================

N = 1000                                # Números de pontos em x (pontos na grid)
xmin, xmax = (-200.,200.)
L = 100.                                # Tamanho do poço - angstrom
L_admensional  = L/A0                   # Tamanho do poço - adimensional
x,dx = np.linspace(xmin/A0,xmax/A0, N, retstep = True, dtype = typem)

print(f'O tamanho do poço é {L:.1f} angstron')
print(f'O passo na Grid é {dx*A0:.1f} angstron')

#===========================================
#>>>>>>>>>>>>>>>>>POTENCIAL<<<<<<<<<<<<<<<<<
#===========================================

V0 = 0.1/Ry                         # Valor do potecial
V = np.zeros(N, dtype = typem)      # Cria um potencial nulo em toda grid

for i in range(N):
    if x[i] <= -L_admensional/2 or x[i] >=L_admensional/2:
        V[i] = V0
    
print(f'O valor do potencial é {V0*Ry:.1f} eV') 

#===========================================
#>>>>>>>>>>>>>>>>HAMILTONIANA<<<<<<<<<<<<<<<
#===========================================

def hamiltoniana(V, N):
    A = (2/dx**2) + V    # Valor da diagonal principal
    B = (-1)/dx**2       # Valor fora da diagonal principal
    H = diags([B, A, B], [-1, 0, 1], shape = (N, N), dtype = typem).toarray()
    return H
H = hamiltoniana(V, N)

#===========================================
#>>>>>>>>>>>>>>>DIAGONALIZAÇÃO<<<<<<<<<<<<<<
#===========================================

def diagonaliza(H, k):
    E, psi = eigsh(H, k=k, which = 'SM', return_eigenvectors = True) # which ='SM' Encontra autovalores e autoestados de 'H' e os retorna em ordem crescente
    return E*Ry, psi

E, psi = diagonaliza(H, k)

#===========================================
#>>>>>>>>>>>>>>>>NORMALIZAÇÃO<<<<<<<<<<<<<<<
#===========================================

def normaliza(psi, k):
    for i in range(0, k):
        integral = np.sum(np.abs(psi[:,i])**2)*dx
        psi[:, i] = psi[:, i]/np.sqrt(integral)
    return psi

norma = normaliza(psi, k)

#===========================================
#>>>>>>>>>DENSIDADE DE PROBABILIDADE<<<<<<<<
#===========================================

prob = np.abs(norma)**2

#===========================================
#>>>>>>GRÁFICOS COM DEN. PROBABILIDADE<<<<<<
#===========================================
fig = plt.figure(figsize=(10,6))
ax1 = fig.add_subplot(111)

for j in range(0, k-1):
    ax1.plot(x*A0, (E[j+1]-E[j])*0.9*prob[:,j]/np.max(prob[:,j])+E[j],
            label = f'$|\psi_{j}(x)|^{2}$')
    
ax1.set_xlabel('x (Angstrom)', fontsize = 12)
ax1.set_ylabel(r'$Energia\; (meV)$', fontsize = 12)
legendh1,labels1=ax1.get_legend_handles_labels()
plt.legend(legendh1,labels1)
plt.savefig("densidade_de_probabilidade.png")
plt.show()

#=============================================
#GRÁFICOS COM DEN. PROBABILIDADE COM POTENCIAL
#=============================================

fig = plt.figure(figsize=(10,6))
ax1 = fig.add_subplot(111)
for j in range(0, k-1):
    ax1.plot(x*A0, (E[j+1]-E[j])*0.9*prob[:,j]/np.max(prob[:,j])+E[j],
            label = f'$|\psi_{j}(x)|^{2}$')
    
ax1.set_xlabel('x (Angstrom)', fontsize = 12)
ax1.set_ylabel(r'$Energia\; (eV)$', fontsize = 12)

plt.plot(x*A0,V*Ry,color="Gray",label="V(x)")
legendh1,labels1=ax1.get_legend_handles_labels()
plt.legend(legendh1,labels1, fontsize = 12)
plt.savefig("densidade_de_probabilidade_com_potencial.png")
plt.show()

def info(H,k):
    print('')
    print('Hamiltonian info')
    print('Número de autoestados =', k)
    print('Formato do Hamiltoniano = ', H.shape)
    print('Tamanho do Hamiltoniano = ', H.size)
    print('')
    
    # Printa valores
    print('System info')
    print('N =',N )
    print('dx ( bohr ) =',dx ,'; dx ( ang ) =', dx*A0 )
    print(f'Potencial = {V0*Ry} (meV)')
    print('L (ang) =',L )
    print('')
    print('Energias')
    for n in range(0, k): # Printa as energias em ordem crescente
        print("E[{}] = {:9.4f} (meV)".format(n,1e3*E[n]))
    return
dados = info(H, k)
print(dados)