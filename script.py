### RHS Right Hand Side ###

# parametros para AC-DC

signal = 0.1
#tasa de produccion basal  
alpha_x = 0.015
alpha_y = 0.3
#señal de activacion (son similares)
beta_x = 6.1
beta_y = 5.7
#tasa de degradacion relativa
#similar degradation rate was observed for all three proteins 
delta_y = 1.05
delta_z = 1.04
#represion
#All the optimizations returned a difference of at least one order 
#of magnitude between the different repression magnitudes (zx < xy < yz < xz)
z_x = 0.00013
x_y = 0.0008
x_z = 0.012
y_z = 0.0012

#exponentes de hill
# no oscillations were found when Hill exponents n = 2 were used, but a small 
#increase of only one of the Hill exponents provided sufficient non-linearity to observe oscillations.
n_zx = 2.3
n_xy = 2
n_xz = 2
n_yz = 2


def RHS(X, P, circuit = 'circuit_name', signal = 0):
    n = X.shape[0] #numero de genes
    G = np.zeros(n) #arreglo de funciones de las interacciones de los genes
        
    if circuit == 'Repressilator':
        for i in range(n):
            j = i+1 if i<n-1 else 0 # el indice del gen que reprime al gen iesimo
            G[i] = P[i][3] * (P[i][0]/(1 + (X[j]/P[i][1])**P[i][2]) - X[i])  

    elif circuit == 'AC-DC' and n == 3:

        G[0] = ((P[0][0] + P[0][1]*signal)/(1+signal+(X[2]/P[2][1])**(P[2][2])) - X[0])
        G[1] = ((P[1][0] + P[1][1]*signal)/(1+signal+(X[0]/P[0][2])**(P[0][4])) - P[1][2]*X[1]) 
        G[2] = (1/((1+(X[0]/P[0][3])**(P[0][5]))+((X[1]/P[1][3])**(P[1][4]))) - P[2][0]*X[2])
    
    
    return G
         
def euler(RHS, X, P, dt):
    return X + RHS(X, P)*dt


## REPRESSILATOR ##
# configuraciones de sistema
# condiciones iniciales de los genes 
X = np.random.random(N) #prueba
N = 3 #numero de genes
#Parametros - provisorio hasta tener el archivo de config
p1 = [10, 9, 3, 2] #parametros gen1 - alfa, beta, gamma, delta 
p2 = [8, 1, 0.5, 6] #parametros gen2 - alfa, beta, gamma, delta
p3 = [6, 5, 0.5, 0.5]

P = [p1, p2, p3]

T = 1000
dt = 0.01

Xt = [X]
for i in range(T):
    Xt.append(euler(RHS, Xt[-1], P ,dt, circuit = 'Repressilator')) 

Xt = np.array(Xt)

## AC-DC ## 
# configuraciones de sistema

X = np.array([0, 0, 0]) #condiciones iniciales
#listas de parametros
p1 = [alpha_x, beta_x, x_y, x_z, n_xy, n_xz]
p2 = [alpha_y, beta_y, delta_y, y_z, n_yz]
p3 = [delta_z, z_x, n_zx]
P = [p1, p2, p3]

T = 1000
dt = 0.1

Xt = [X]
for i in range(T):
    Xt.append(euler(RHS, Xt[-1], P ,dt, circuit = 'AC-DC', signal = 0.1)) 

Xt = np.array(Xt)
