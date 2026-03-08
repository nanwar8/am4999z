import numpy as np
import matplotlib.pyplot as plt

a = 1      
k1 = 12.5 
G1 = (72.8)**2   
G2 = 10  
d = 0.09

b = 0      
k2 = 10  
g = 10   

j_vals = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

# reduced parameters
G3 = (G2*b)/g
G4 = (G2*k2)/g
G5 = G1 + G3

xs = np.linspace(0,100,1000)

for j in j_vals:

    ys = (k1*(xs**2))/(G5 + xs**2 + G4*xs) - d*xs + j*(a/(d+j))

    plt.figure()
    plt.plot(xs, ys)
    plt.axhline(0, linestyle="--")
    
    plt.ylim(-1,1)
    plt.xlim(0,100)

    plt.title(f"j = {j}")
    plt.xlabel("p2")
    plt.ylabel("f(p2)")

    plt.show()