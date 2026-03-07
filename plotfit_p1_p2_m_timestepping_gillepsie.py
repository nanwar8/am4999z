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

# calculating the individual rates
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
    

# updating the number of p1, p2, and m
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

# implement gillepsie algorithm, put in initial variables and constants from case 3
# produce list of times, p1, p2, and m concentrations by computing all reaction rates; next rate depends on sampling of tau
# meaning times are NOT uniform, times are dependent on the spacing between saved times steps in EVERY time step
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
   

def gillepsie_uniformtimes(p10, p20, m0, a, k1, G1, G2, d, b, k2, g, j, t_max, t_step):

    p1_list, p2_list, m_list, t_list = gillepsie_two_boxes(p10, p20, m0, a, k1, G1, G2, d, b, k2, g, j, t_max)
    
    t_list=np.asarray(t_list)
    p1_list=np.asarray(p1_list)
    p2_list=np.asarray(p2_list)
    m_list=np.asarray(m_list)
    
    t_unigrid = np.arange(0,t_max+t_step,t_step)
    
    uniform_idx = np.searchsorted(t_list, t_unigrid, side="right") -1
    uniform_idx=np.clip(uniform_idx,0, len(t_list)-1)
    
    p1_unigrid = p1_list[uniform_idx]
    p2_unigrid = p2_list[uniform_idx]
    m_unigrid = m_list[uniform_idx]
    return p1_unigrid.tolist(), p2_unigrid.tolist(), m_unigrid.tolist(), t_unigrid.tolist()
   
def gillespie_two_boxes_on_grid(p10, p20, m0, a, k1, G1, G2, d, b, k2, g, j, t_max, t_step):
    t = 0.0
    p1, p2, m = p10, p20, m0

    t_grid = np.arange(0, t_max + t_step, t_step)
    p1_grid = np.empty_like(t_grid, dtype=int)
    p2_grid = np.empty_like(t_grid, dtype=int)
    m_grid  = np.empty_like(t_grid, dtype=int)

    grid_i = 0
    p1_grid[grid_i], p2_grid[grid_i], m_grid[grid_i] = p1, p2, m

    # precompute next sample time
    next_sample_t = t_grid[grid_i + 1] if grid_i + 1 < len(t_grid) else None

    while t < t_max:
        rates = rate_comp_p_m(p1, p2, m, a, k1, G1, G2, d, b, k2, g, j)
        lam = sum(rates)
        if lam <= 0:
            # fill remaining grid with last value
            while grid_i + 1 < len(t_grid):
                grid_i += 1
                p1_grid[grid_i], p2_grid[grid_i], m_grid[grid_i] = p1, p2, m
            break

        tau = -math.log(np.random.rand()) / lam
        t_next = t + tau

        # if we crossed one or more sample times, write out the last state
        while next_sample_t is not None and t_next >= next_sample_t:
            grid_i += 1
            p1_grid[grid_i], p2_grid[grid_i], m_grid[grid_i] = p1, p2, m
            next_sample_t = t_grid[grid_i + 1] if grid_i + 1 < len(t_grid) else None
            if next_sample_t is None:
                break

        if t_next > t_max:
            break

        # advance time and apply one event
        t = t_next
        event = choosing_events(rates)
        p1, p2, m = update_numbers(event, p1, p2, m)

    return p1_grid.tolist(), p2_grid.tolist(), m_grid.tolist(), t_grid.tolist()
    
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

t_max = 100000
t_step = 20
j=0.3

p1_list0, p2_list0, m_list0, t_list0 = gillepsie_uniformtimes(p10, p20, m0, a, k1, G1, G2, d, b, k2, g, j, t_max, t_step)
p1_list, p2_list, m_list, t_list = gillespie_two_boxes_on_grid(p10, p20, m0, a, k1, G1, G2, d, b, k2, g, j, t_max, t_step)

# j_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
# j_vals = [0.6, 0.7, 0.8, 0.9, 1.0]

t_max = 100000
timestep = 20

binsize_p1 = 15
binsize_p2 = 35
binsize_m = 45

for j in j_vals:
    p1_list, p2_list, m_list, times = gillespie_two_boxes_on_grid(
        p10, p20, m0, a, k1, G1, G2, d, b, k2, g, j, t_max, timestep
    )


    t_plotstart = 10000 
    t_axis = np.array(times)  # already uniform
    start_idx = int(t_plotstart / timestep)

# p1
    plt.figure()
    plt.plot(t_axis[start_idx:], p1_list[start_idx:])
    plt.title(f"time series of $\\bar{{p1}}$ (Δt = {timestep}, j = {j})")
    plt.xlabel("time")
    plt.ylabel("number of p1 molecules")

    plt.figure()
    plt.hist(p1_list[start_idx:], bins=binsize_p1)
    plt.title(f"histogram of $\\bar{{p1}}$ (Δt = {timestep}, j = {j}, binsize={binsize_p1})")
    plt.xlabel("$\\bar{p1}$ per time bin")
    plt.ylabel("count")

# p2
    plt.figure()
    plt.plot(t_axis[start_idx:], p2_list[start_idx:])
    plt.title(f"time series of $\\bar{{p2}}$ (Δt = {timestep}, j = {j})")
    plt.xlabel("time")
    plt.ylabel("number of p2 molecules")

    plt.figure()
    plt.hist(p2_list[start_idx:], bins=binsize_p2)
    plt.title(f"histogram of $\\bar{{p2}}$ (Δt = {timestep}, j = {j}, binsize={binsize_p2})")
    plt.xlabel("$\\bar{p2}$ per time bin")
    plt.ylabel("count")


# m
    plt.figure()
    plt.plot(t_axis[start_idx:], m_list[start_idx:])
    plt.title(f"time series of $\\bar{{m}}$ (Δt = {timestep}, j = {j})")
    plt.xlabel("time")
    plt.ylabel("number of m molecules")
    
    plt.figure()
    plt.hist(m_list[start_idx:], bins=binsize_m)
    plt.title(f"histogram of $\\bar{{m}}$ (Δt = {timestep}, j = {j}, binsize={binsize_m})")
    plt.xlabel("$\\bar{m}$ per time bin")
    plt.ylabel("count")


    plt.show()