### RHS Right Hand Side ###

# configuraciones de sistema
N = 4 #numero de genes
#Parametros - provisorio hasta tener el archivo de config
p1 = [10, 2, 0.7, 2] #parametros gen1 - alfa, beta, gamma, delta 
p2 = [13, 1, 0.3, 0.5] #parametros gen2 - alfa, beta, gamma, delta
p3 = [13, 1, 0.3, 0.5]
p4 = [13, 1, 0.3, 0.5]
P = [p1, p2, p3, p4]


# condiciones iniciales de los genes 
X = np.random.random(N) #prueba

#repressilator
def RHS(X, P):
    n = X.shape[0] #numero de genes
    G = np.zeros(n) #arreglo de funciones de las interacciones de los genes
    for i in range(n):
        j = i+1 if i<n-1 else 0 # el indice del gen que reprime al gen iesimo
        G[i] = P[i][3] * (P[i][0]/(1 + (X[j]/P[i][1])**P[i][2]) - X[i])  
    
    return G
         
def euler(RHS, X, P, dt):
    return X + RHS(X, P)*dt


T = 1000
dt = 0.01

Xt = [X]
for i in range(T):
    Xt.append(euler(RHS, Xt[-1], P ,dt)) 

Xt = np.array(Xt)

plt.plot(Xt)