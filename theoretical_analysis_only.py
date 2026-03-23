import numpy as np
import matplotlib.pyplot as plt

a = 1      
k1 = 12.5 
G1 = (72.8)**(2)   
G2 = 10  
d = 0.09

b = 0      
k2 = 10  
g = 10   

# initial conditions
p10 = 100 
p20 = 50    
m0 = 30

j_vals = [0.1,0.4]

# reduced parameters
G3 = (G2*b)/g
G4 = (G2*k2)/g
G5 = G1 + G3

t_max = 200000
timestep = 1

for j in j_vals:
    # theoretical
    xs = np.linspace(0,200,1000)
    ys = (k1*(xs**2))/(G5 + xs**2 + G4*xs) - d*xs + j*(a/(d+j))

    plt.figure()
    plt.plot(xs,ys)
    plt.xlim(0,140)
    plt.ylim(-2,2)
    plt.title(f"j = {j}")
    plt.ylabel('f(p2)')
    plt.axhline(0)

    plt.show()