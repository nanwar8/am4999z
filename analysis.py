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

j_vals = [0.1]

# reduced parameters
G3 = (G2*b)/g
G4 = (G2*k2)/g
G5 = G1 + G3

t_max = 200000
timestep = 1

for j in j_vals:
    p1_list, p2_list, m_list, times = gillepsie_two_boxes(p10, p20, m0, a, k1, G1, G2, d, b, k2, g, j, t_max)

    #averages_p1, counts_p1 = timestepavg(times, p1_list, timestep)
    averages_p2, counts_p2 = timestepavg(times, p2_list, timestep)
    #averages_m, counts_m = timestepavg(times, m_list, timestep)

    p2_arr=np.array(averages_p2)
    
    p2result = []

    for i in range(len(p2_arr)-1):
        current = p2_arr[i]
        diff = p2_arr[i+1] - p2_arr[i]
        p2result.append([current, diff])
        
    p2finresult = np.array(p2result)
    p2finresult = p2finresult[p2finresult[:,0].argsort()] #sorted in ascending order from 1st col
    
    interval=20
    
    x = p2finresult[:,0]   # 0th column
    y = p2finresult[:,1]   # 1st column
    
    result=[]
    
    start=np.min(x)
    end = np.max(x)
    
    for left in np.arange(start, end + interval, interval):
        right = left + interval 
        mask = (x >= left) & (x < right)
        
        if np.any(mask):
            x_center = (left + right) / 2
            y_avg = np.mean(y[mask])
            result.append([x_center, y_avg])
            
    result = np.array(result)
    
    x=result[:,0]
    y=result[:,1]
    # theoretical
    xs = np.linspace(0,100,1000)
    ys = (k1*(xs**2))/(G5 + xs**2 + G4*xs) - d*xs + j*(a/(d+j))

    plt.figure()
    plt.plot(x,y)
    plt.plot(xs,ys)
    plt.ylim(np.min(y),np.max(y))
    plt.xlim(0,80)
    plt.title(f"j = {j}")
    plt.axhline(0)

    plt.show()