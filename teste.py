import numpy as np 
from scipy import constants as cte

def potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N):
	Ry = cte.value(u'Rydberg constant times hc in eV')
	A0 = cte.value(u'Bohr radius')/cte.value(u'Angstrom star')
	V00 = VPOCO
	V0 = (V00*1e-3)/Ry
	V_POROS = VPOROS/Ry
	LX = Largura
	LPX = LX
	RMAX = 5
	LPY = Largura
	Largura_adm  = Largura/A0
	LRX = (LPX - 2*RMAX)/A0
	LRY = (LPY - 2*RMAX)/A0
	V = np.zeros((N,N), dtype = np.float64)
	VTESTE = V
	
	if poros == 1:
		for i in range(0, num_poros):
			NUMTESP = 0
			while NUMTESP == 0:
				NSUPERPONDO = 0
				VTESTE = VTESTE*0
				XP = np.random.rand()
				XP = ((XP-0.5))*LRX
				YP = np.random.rand()
				YP = ((YP-0.5))*LRY
				
				for IX in range(len(x)):
					for IY in range(len(y)):
						RP = np.sqrt((x[IX]-YP)**2 + (y[IY]-XP)**2)
						if(RP<RMAX):
							VTESTE[IX,IY] = V_POROS
						if(V[IX,IY] + VTESTE[IX,IY] > V_POROS):
							NSUPERPONDO = NSUPERPONDO + 1
						if(y[IY] <= -Largura_adm/2 or y[IY] >= Largura_adm/2):
							V[IX,IY] = V0
						if(x[IX] <= -Largura_adm/2 or x[IX] >= Largura_adm/2):
							V[IX,IY] = V0			
				if(NSUPERPONDO == 0):
					V = V + VTESTE
					NUMTESP = 1
		potencial = V
		
	elif poros == 0:
		for IX in range(len(x)):
			for IY in range(len(y)):
				if(x[IX] <= -Largura_adm/2 or x[IX] >= Largura_adm/2):
					V[IX,IY] = V0
				if(y[IY] <= -Largura_adm/2 or y[IY] >= Largura_adm/2):
					V[IX,IY] = V0
		potencial = V

	return potencial

#=================================================
#=========== Função que descreve a massa =========
#=================================================
import numpy as np 
from scipy import constants as cte

def massa(x,y,Largura,V, VPOROS, VPOCO):
  A0 = cte.value(u'Bohr radius')/cte.value(u'Angstrom star')
  Ry = cte.value(u'Rydberg constant times hc in eV')
  V00 = VPOCO
  V0 = (V00*1e-3)/Ry
  V_POROS = VPOROS/Ry
  m1 = 0.041  # InGaAs
  m2 = 0.075      # InAlAS
  m3 = 1      # Vácuo
  L = Largura
  Largura_adm = L/A0

  m_eff = m1*np.ones((len(x),len(y)), dtype = np.float64)    #effective mass
  
  for i in range(len(x)):
    for j in range(len(y)):
      #if V[i,j] == V0:
      #  m_eff[i,j] = m2
      if V[i,j] == V_POROS: 
        m_eff[i,j] = m3

# teste mostrado p/ teldo: massa igual no poço, potencial infinito
 # for i in range(len(x)):
 #   for j in range(len(y)):
 #     if  x[i] <= -Largura_adm/2 or x[i] >= Largura_adm/2:
 #       m_eff[i,j] = m2
 #     if  y[j] <= -Largura_adm/2 or y[j]>= Largura_adm/2:
  #      m_eff[i,j] = m2

  return m_eff

import numpy as np                           # Biblioteca de matrizes
import matplotlib.pyplot as plt              # Biblioteca de gráficos
from scipy.sparse.linalg import eigsh        # Biblioteca para Diagonalizar a matriz hamiltoniana                   # Biblioteca de performace
from scipy.sparse import diags               # Biblioteca de cria uma matriz tridiagonal (sparse)
import matplotlib.pyplot as plt
from scipy import constants as cte
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d

def hamiltoniana(V, mu, dx, dy, N):
  #mu --------> mu = m0/m*(x,y) ----> m*(x,y)=0.041m0
  #poteff ---------> Potencial adimensional
  #dx,dy --------> Passo - adimensional
  #N ---------> Número de pontos

    A = np.zeros((N,N), dtype = np.float64)        #Coeficiente da diagonal principal
    # A[i,j] = ((mu[i+1,j]+2*mu[i,j]+mu[i-1,j])/(2*(dx**2))) + (((mu[i,j+1]+2*mu[i,j]+mu[i,j-1]))/(2*(dy**2))) + potef[i,j]

    for i in range(N):
        for j in range(N):
            if (i == 0  and j == 0): #3
                A[0,0] = ((mu[1,0]+2*mu[0,0]+mu[0,0])*(1/(2*(dx**2)))) + ((mu[0,1]+2*mu[0,0]+mu[0,0])*(1/(2*(dy**2)))) + V[0,0]
            if (i == N-1 and j == N-1): #3
                A[N-1,N-1]= ((mu[0,N-1]+2*mu[N-1,N-1]+mu[N-2,N-1])*(1/(2*(dx**2)))) + ((mu[N-1,0]+2*mu[N-1,N-1]+mu[N-1,N-2]))*(1/(2*(dy**2))) + V[N-1,N-1]
#========================== PEGA TODOS OS TERMOS DE i=0 e todos os valores de j===================================================
            if(i == 0 and j!= N-1):   #3
                A[0,j] = ((mu[1,j]+2*mu[0,j]+mu[0,j])*(1/(2*(dx**2)))) + ((mu[0,j+1]+2*mu[0,j]+mu[0,j-1]))*(1/(2*(dy**2))) + V[0,j]

            if(i != 0 and j == N-1):  #3
                A[i,N-1] = ((mu[i,N-1]+2*mu[i,N-1]+mu[i-1,N-1])*(1/(2*(dx**2)))) + ((mu[i,0]+2*mu[i,N-1]+mu[i,N-2])*(1/(2*(dy**2)))) + V[i,N-1]####

            if(i != 0 and j == 0):
                A[i,0] = ((mu[0,0]+2*mu[0,0]+mu[i-1,0])*(1/(2*(dx**2)))) + ((mu[i,1]+2*mu[i,0]+mu[i,0])*(1/(2*(dy**2)))) + V[i,0]
               
            if (i != N-1  and j !=N-1):
                A[i,j] = ((mu[i+1,j]+2*mu[i,j]+mu[i-1,j])*(1/(2*(dx**2)))) + ((mu[i,j+1]+2*mu[i,j]+mu[i,j-1])*(1/(2*(dy**2)))) + V[i,j]
               
            if(i == N-1 and  j != 0):
                A[N-1,j] = ((mu[0,0]+2*mu[N-1,0]+mu[N-2,0])*(1/(2*(dx**2)))) + ((mu[N-1,0]+2*mu[N-1,0]+mu[N-1,j-1])*(1/(2*(dy**2)))) + V[N-1,j]

            if(i == N-1 and  j == 0):#3
                A[N-1,0] = ((mu[0,0]+2*mu[N-1,0]+mu[N-2,0])*(1/(2*(dx**2)))) + ((mu[N-1,1]+2*mu[N-1,0]+mu[N-1,0])*(1/(2*(dy**2)))) + V[N-1,0]
            
            if(i!=0 and j!=0):
                A[i,j] = ((mu[0,j]+2*mu[i,j]+mu[i-1,j])*(1/(2*(dx**2)))) + ((mu[i,0]+2*mu[i,j]+mu[i,j-1])*(1/(2*(dy**2)))) + V[i,j]
            
            if(i==0 and j==N-1):
                A[0,N-1] = ((mu[1,N-1]+2*mu[0,N-1]+mu[0,N-1])*(1/(2*(dx**2)))) + ((mu[0,0]+2*mu[0,N-1]+mu[0,N-2])*(1/(2*(dy**2)))) + V[0,N-1]
            
    a = A.flatten()
#===================================================================
    B = np.zeros((N,N), dtype = np.float64)
    # B[i,j] = -(mu[i,j]+mu[i+1,j])*(1/(2*(dx**2)))
    
    for i in range(N):
        for j in range(N):
            if (i == 0 and j == 0): #3
                B[0,0] = -(mu[0,0]+mu[1,0])*(1/(2*(dx**2)))
            if(i == N-1 and j == N-1 ): 
                B[N-1,N-1] = -(mu[N-1,N-1]+mu[0,N-1])*(1/(2*(dx**2)))
                
            if (i == 0  and j !=N-1):  #3
                B[0,j] = -(mu[0,j]+mu[1,j])*(1/(2*(dx**2)))
            
            if(i == N-1 and  j != 0):
                B[N-1,j] = -(mu[N-1,j]+mu[0,j])*(1/(2*(dx**2)))  # Zerou o índice i+1
    #--------------------------------------------------------------------------            
            if (i != 0  and j ==0):
                B[i,0] = -(mu[i,0]+mu[0,0])*(1/(2*(dx**2)))  # Zerou o índice i+1
    #--------------------------------------------------
            if (i != N-1  and j !=N-1):
                B[i,j] = -(mu[i,j]+mu[i+1,j])*(1/(2*(dx**2)))    
    b =B.flatten()
#===================================================================

    B1 = np.zeros((N,N), dtype = np.float64)
    # B[i,j] = -(mu[i,j]-mu[i-1,j])*(1/(2*(dx**2)))

    for i in range(N):
        for j in range(N):
            if (i == 0 and j == 0): #3
                B1[0,0] = -(mu[0,0]+mu[0,0])*(1/(2*(dx**2)))
            if(i == N-1 and j == N-1 ): 
                B1[N-1,N-1] = -(mu[N-1,N-1]+mu[N-2,N-1])*(1/(2*(dx**2)))
                
            if (i == 0  and j !=N-1):  #3
                B1[0,j] = -(mu[0,j]+mu[0,j])*(1/(2*(dx**2)))
            
            if(i == N-1 and  j != 0):
                B1[N-1,j] = -(mu[N-1,j]+mu[N-2,j])*(1/(2*(dx**2)))  # Zerou o índice i+1
    #--------------------------------------------------------------------------            
            if (i != 0  and j ==0):
                B1[i,0] = -(mu[i,0]+mu[i-1,0])*(1/(2*(dx**2)))  # Zerou o índice i+1
    #--------------------------------------------------
            if (i != N-1  and j !=N-1):
                B1[i,j] = -(mu[i,j]+mu[i-1,j])*(1/(2*(dx**2)))
           
    b1 =B.flatten()
    #===================================================================

    C = np.zeros((N,N), dtype = np.float64)
    #C[i,j] = -(mu[i,j]+mu[i,j+1])*(1/(2*(dy**2)))

    for i in range(N):
        for j in range(N):
            if (i==0 and j == 0):  #3
                C[0,0] = -(mu[0,0]+mu[0,1])*(1/(2*(dy**2)))
            if(i == N-1 and j == N-1 ):   #3
                C[N-1,N-1] = -(mu[N-1,N-1]+mu[N-1,0])*(1/(2*(dy**2)))
                
            if (i == 0  and j == N-1): 
                C[0,N-1] = -(mu[0,N-1]+mu[0,0])*(1/(2*(dy**2)))
                    
            if(i !=0 and  j == N-1):   #3
                C[i,N-1] = -(mu[i,N-1]+mu[i,0])*(1/(2*(dy**2)))
                
            if (i != N-1  and j !=N-1):
                C[i,j] = -(mu[i,j]+mu[i,j+1])*(1/(2*(dy**2)))
    
    c = C.flatten()

        #===================================================================
    C1 = np.zeros((N,N), dtype = np.float64)
    #C[i,j] = -(mu[i,j]+mu[i,j-1])*(1/(2*(dy**2)))

    for i in range(N):
        for j in range(N):
            if (i==0 and j == 0):  #3
                C1[0,0] = -(mu[0,0]+mu[0,0])*(1/(2*(dy**2)))
            if(i == N-1 and j == N-1 ):   #3
                C1[N-1,N-1] = -(mu[N-1,N-1]+mu[N-1,N-2])*(1/(2*(dy**2)))
                
            if (i == 0  and j == N-1): 
                C1[0,N-1] = -(mu[0,N-1]+mu[0,N-2])*(1/(2*(dy**2)))
                    
            if(i !=0 and  j == N-1):   #3
                C1[i,N-1] = -(mu[i,N-1]+mu[i,N-2])*(1/(2*(dy**2)))
                
            if (i != N-1  and j !=N-1):
                C1[i,j] = -(mu[i,j]+mu[i,j-1])*(1/(2*(dy**2)))
    
    c1 = C.flatten()
    #===================================================================
    H = diags([c,b,a,b,c], [-N, -1, 0, 1, N], shape =(N**2, N**2), dtype = np.float64).tocsr()
    
    return  H

    from scipy import constants as cte

def diagonaliza(H, k, dx, dy):
    Ry = cte.value(u'Rydberg constant times hc in eV')
    print (' ')
    print ('Hamiltonian Info ')
    print ('Number of states / eigenvalues =',k )
    print ('Matrix shape =', H.shape )
    print ('Matrix n. of elements =', H.size )
    print (' ')
    En, psi = eigsh(H, k=k, which = 'SM', return_eigenvectors = True) # which ='SM' Encontra autovalores e autoestados de 'H' e os retorna em ordem crescente

    for i in range(0, k):
        integral = np.sum(np.abs(psi[:,i])**2)*dx*dy
        psi[:, i] = psi[:, i]/np.sqrt(integral)
    return En*Ry*1e3, psi

def config(V,ex,xmin,xmax):
    fig = plt.figure(figsize=(10,8)) #terrain_r, YlGnBu
    ax1 = fig.add_subplot(111)
    plt.imshow(V*Ry*1e3, extent = ex, origin='lower', cmap = 'YlGnBu', interpolation = 'nearest')    
    ax1.set_xlabel('x (Angstrom)', fontsize = 14)
    ax1.set_ylabel('y (Angstrom)', fontsize = 14)
    plt.yticks(range(xmin,xmax+50,50))
    plt.xticks(range(xmin,xmax+50,50))
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    cb = plt.colorbar()
    cb.set_label(r'V(x,y)$\; (me$V)', fontsize = 14)
    return


def plot_grid(configuracao, Largura, c):
    plt.grid( color='0.95',alpha=0.05)
    plt.savefig(f"Contagem - {c} - potencial-L={Largura}-com-grid.png")
    plt.savefig(f"Contagem - {c} - potencial-L={Largura}-com-grid-.svgz")
    plt.show()
    return

   
def plot_sem_grid(configuracao, Largura, c):
    plt.savefig(f"Contagem - {c} - potencial-L={Largura}.png")
    plt.savefig(f"Contagem - {c} - potencial-L={Largura}.svgz")
    plt.show()
    return 

from scipy import constants as cte

def alfa(alpha,x,y,V, VPOROS, VPOCO):
    Ry = cte.value(u'Rydberg constant times hc in eV')
    V00 = VPOCO
    V0 = (V00*1e-3)/Ry
    V_POROS = VPOROS/Ry 
    for i in range(len(x)):
        for j in range(len(y)):
            if V[i,j] == 0:
                alpha[i,j] = 0.2
            if V[i,j] == V_POROS: 
                alpha[i,j] = 0
            if V[i,j] == V0: 
                alpha[i,j] = 0

    return alpha

from scipy import constants as cte

def dados_funçao_de_onda(alf4, V, psi, N, k, x,X,Y, Largura):
    A0 = cte.value(u'Bohr radius')/cte.value(u'Angstrom star') 
    Ry = cte.value(u'Rydberg constant times hc in eV')

    x1 = X*A0
    y1 = Y*A0
    V1 = V*Ry*1e3
    V2 = V*Ry*alf4
    
    for i in range(k):
        wave = psi.T[i].reshape((N,N))
        wave2D = wave**2
        wave2D_max = np.amax(wave2D)
        funcao_de_onda = wave2D/wave2D_max

        for j in range(len(x)):
            potencial = V1.flatten()
            potencial_alfa = V2.flatten()
            x_1 = x1.flatten()
            y_1 = y1.flatten()
            fun = funcao_de_onda.flatten()
            data = np.column_stack([x_1, y_1, potencial_alfa, potencial, fun])
            np.savetxt(f'dados_L={Largura}_psi{i}_com_potencial.dat', data)
    return

def plot_psi(k,V,ex,xmin,xmax,psi,alf4,Largura,i, c):
    for z in range(len(c)):
        for n in range(len(i)):
            for j in range(k):
                wave = psi.T[j].reshape((N,N))
                wavef2D = wave**2
                wavef2D_max = np.amax(wavef2D)#viridis, YlGnBu,terrain_r
                
                fig = plt.figure(figsize=(10,8))
                ax1 = fig.add_subplot(111)
                plt.imshow(wavef2D/wavef2D_max , extent = ex ,cmap =  'jet', origin='lower', interpolation = 'nearest')
                cb = plt.colorbar()
                cb.set_label(label = r'$ \frac{|\psi_{n}(x,y)|^{2}}{|\psi_{max}(x,y)|^{2}}$', fontsize = 14)
                plt.imshow(V*Ry*1e3, extent = ex, alpha= alf4, origin ='lower', cmap = 'YlGnBu', interpolation = 'nearest')
                ax1.set_xlabel('x (Angstrom)', fontsize = 14)
                ax1.set_ylabel('y (Angstrom)', fontsize = 14)
                plt.yticks(range(xmin,xmax+50,50))
                plt.xticks(range(xmin,xmax+50,50))
                plt.xticks(fontsize = 14)
                plt.yticks(fontsize = 14)
                plt.savefig(f"Contagem - {c} - L={Largura} - snap-shot {j}.png")
                plt.savefig(f"Contagem - {c} - L={Largura} - snap-shot {j}.svgz")
                plt.show()

def dados_plot(Largura, En, c):
        nome = f'dados - Contagem - {c} - LxE_L={Largura}.dat'
        file = open(nome,'w')
        file.write('Lx=Ly       E0           E1           E2           E3           E4           E5\n')
        file.write(f'{Largura:5.5f}    {En[0]:5.5f}    {En[1]:5.5f}    {En[2]:5.5f}    {En[3]:5.5f}    {En[4]:5.5f}    {En[5]:5.5f}\n')
        file.close()   

import numpy as np                           # Biblioteca de matrizes
from scipy import constants as cte
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
A0 = cte.value(u'Bohr radius')/cte.value(u'Angstrom star')  # Raio de Bohr ---> (Angstroms)
Ry = cte.value(u'Rydberg constant times hc in eV')

N = 1000  
xmin, xmax = (-150,150)      
ymin, ymax = (xmin, xmax)           # Tamanho do poço - adimensional
ex = [xmin, xmax, ymin, ymax]
x, dx = np.linspace(xmin/A0,xmax/A0, N, retstep = True, dtype = np.float64 )  # Define o valor de x e dx
y, dy = (x, dx)
X, Y = np.meshgrid(x,y)
VPOCO = 600
VPOROS = 10
k = 6

alpha = np.zeros((N,N), dtype=np.float64)
contagem = range(8,10)
largura_poço = range(40,210,10)

poros = 1 

for c in contagem: 
    for i in largura_poço:
        if i == 40:
            num_poros = 5
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 50:
            num_poros = 8
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
            
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)

        if i == 60:
            num_poros = 11
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 70:
            num_poros = 15
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 80:
            num_poros = 20
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 90:
            num_poros = 25
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 100:
            num_poros = 32
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 110:
            num_poros = 38
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 120:
            num_poros = 46
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 130:
            num_poros = 53
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 140:
            num_poros = 62
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 150:
            num_poros = 71
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)

        if i == 160:
            num_poros = 81
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 170:
            num_poros = 92
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 180:
            num_poros = 103
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 190:
            num_poros = 115
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)
        
        if i == 200:
            num_poros = 127
            Largura = i

            V = potencial(x, y, VPOROS, VPOCO, Largura, poros, num_poros, N)
                
            massa_eff = massa(x,y,Largura,V, VPOROS, VPOCO)
            mu = np.divide(1,massa_eff)

            configuracao = config(V,ex,xmin,xmax)
            figura_sem_gird = plot_sem_grid(configuracao, Largura, [c])
                    
            H = hamiltoniana(V, mu, dx, dy, N)
            En, psi = diagonaliza(H, k, dx, dy)

            alf4 = alfa(alpha,x,y,V, VPOROS, VPOCO) 
            snap_psi = plot_psi(k,V,ex,xmin,xmax,psi,alf4, Largura, [i], [c])
            
            dados = dados_plot(Largura,En, [c])

            for n in range(0,k):
                print(f'E{n} = {En[n]:.4f}')
            print(i)