import numpy as np
import matplotlib.pyplot as plt
import math
np.random.seed(60)


def exponential_sampling(lam):
    u = np.random.rand()
    tau = -np.log(u) / lam
    return tau

def choosing_events(rates):
    sum_rates = sum(rates)
    rate_probabilities = [r / sum_rates for r in rates] 
    u2 = np.random.rand() 
    
    cumul = 0
    for i, p in enumerate(rate_probabilities): 
        cumul += p 
        if u2 < cumul :
            return i 

def rate_comp_p_m(p1,p2,m,a,k1,G1,G2,d,b,k2,g,j):
    rates=[]
    
    # first box contributors
    r0 = a # add to p1
    rates.append(r0)
    
    r1 = d*p1 # kill p1
    rates.append(r1)
    
    r2 = (k1 * p2**2) / (G1 + p2**2 + G2*m) # add to p2
    rates.append(r2)
    
    # second box contributors
    
    r3 = b # m production
    rates.append(r3)
    
    r4 = k2 * p2 # m production
    rates.append(r4) 
    
    r5 = g * m # kill m
    rates.append(r5) 
    
    r6 = d * p2 # kill p2
    rates.append(r6)

    # diffusion of p1 to p2 
    r7 = j * p1 # turn p1 into p2
    rates.append(r7)
    
    return rates
       
def update_numbers(event, p1, p2, m):
    if event == 0:
        p1 += 1
    elif event == 1:
        if p1 > 0: p1 -= 1
    elif event == 2:
        p2 += 1
    elif event == 3:
        m += 1
    elif event == 4:
        m += 1
    elif event == 5:
        if m > 0: m -= 1
    elif event == 6:
        if p2 > 0: p2 -= 1
    elif event == 7:
        if p1 > 0:
            p1 -= 1
            p2 += 1
    return p1, p2, m

def gillepsie_two_boxes(p10, p20, m0, a, k1, G1, G2, d, b, k2, g, j, t_max):
    t = 0
    p1, p2, m = p10, p20, m0
    
    times = [t]
    p1_list = [p1]
    p2_list = [p2]
    m_list=[m]
    
    event_counts = np.zeros(8) 
    while t < t_max:
      rate_array = rate_comp_p_m(p1,p2,m,a,k1,G1,G2,d,b,k2,g,j)
      lam = sum(rate_array)
      if lam <= 0: 
          break 
      
        
      tau = exponential_sampling(lam)
      t_next = t + tau
      if t_next > t_max: 
          break 
      t = t_next
      
      event = choosing_events(rate_array)
      event_counts[event] += 1
      
      p1, p2,m = update_numbers(event, p1, p2, m)
      
      p1_list.append(p1)
      p2_list.append(p2)
      m_list.append(m)
      times.append(t) 
        ## density check
    return p1_list, p2_list, m_list, times
   
def timestepavg(times, x_list, timestep):
    max_time = max(times)
    binsize = int(math.floor(max_time) / timestep) + 1

    sums = [0] * binsize
    counts = [0] * binsize

    for i in range(len(times)):
        b = int(math.floor((times[i]) / timestep))
        sums[b] += + x_list[i]
        counts[b] += 1
        
    averages = [0] * binsize
    for b in range(binsize):
        if counts[b] != 0:
            averages[b] = sums[b] / counts[b]
    averages[0] = x_list[0]

    return averages, counts

def make_bins(data, bin_width=1, start=0):
    data = np.asarray(data)
    xmax = np.max(data)
    return np.arange(start, xmax + bin_width, bin_width)
### testing ### 
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

j_vals = [0.1, 0.4]

# reduced parameters
G3 = (G2*b)/g
G4 = (G2*k2)/g
G5 = G1 + G3


t_max = 100000
timestep = 1

#binsize = 45

for j in j_vals:
    p1_list, p2_list, m_list, times = gillepsie_two_boxes(p10, p20, m0, a, k1, G1, G2, d, b, k2, g, j, t_max)

    averages_p2, counts_p2 = timestepavg(times, p2_list, timestep)

    t_plotstart = 1000 
    start_idx = int(t_plotstart / timestep) 
    t_axis = np.arange(len(averages_p2)) * timestep

    data_p2 = np.array(averages_p2[start_idx:])
    
# time series 
    plt.figure()
    plt.plot(t_axis[start_idx:], averages_p2[start_idx:])
    plt.title(f"Time Series of p2 (j = {j})")
    plt.xlabel("time")
    plt.ylabel("number of p2 molecules")

# histogram with theoretical curves
    fig, ax = plt.subplots(2, 1, figsize=(6, 8))

# histogram
    ax[0].hist(data_p2, bins=make_bins(data_p2, bin_width=1, start=0))
    ax[0].set_title(f"Histogram of p2 (j = {j})")
    ax[0].set_xlim(0, 140)
    ax[0].set_ylabel("frequency")

# theoretical
    xs = np.linspace(0, 200, 1000)
    ys = (k1*(xs**2))/(G5 + xs**2 + G4*xs) - d*xs + j*(a/(d+j))

    ax[1].plot(xs, ys)
    ax[1].set_xlim(0, 140)
    ax[1].set_ylim(-2, 2)
    ax[1].set_title(f"f(p2), j = {j}")
    ax[1].set_xlabel("number of p2 molecules (p2)")
    ax[1].set_ylabel('dp/dt')
    ax[1].axhline(0)

    plt.tight_layout()
    plt.show()
