# Libraries 
import numpy as np
from scipy.interpolate import interp1d





# Function to size solar collector
def dynamic_power_dependence(SWP,P_in,I_in,t,dt=5,print_results=1,V_feed=[],V_ratio=[],A=[],interp_kind='cubic'):
    """
    Estimates performance of desalination system under time-varying solar intensity using
    stationary simulation results. 

    Sizes optimal solar collector footprints in terms of efficiency and production given:
    1. Optimized efficiency "SWP" vs. power "P_in" curve of a given system
    2. Time-varying solar intensity profile "I_in" vs. "t"
    3. Time discretization "dt"

    Inputs:
    - SWP (L/kWh): array containing maximum SWP below boiling of a given system vs. powers in P_in
    - P_in (W): powers at which the optimal SWP's were evaluated
    - I_in (W/m^2): solar intensities vs. times in t
    - t (h): times at which the intensities were measured

    Optional Inputs:
    - dt (min): time discretization/interval over which to evaluate average production rate (DEFAULT 5 min)
    - print_results: print aggregate results as text if true (DEFAULT True)
    - V_feed (L/h): array containing optimized feed flow rates vs. powers in P_in (DEFAULT EMPTY)
    - V_ratio (1): array containing optimized flow rate ratios vs. powers in P_in (DEFAULT EMPTY)

    ***dt MUST be a multiple of 5 min***

    Outputs:
    - A (m^2): array containing all tested collector sizes
    - W (L): array containing all total production values vs. collector size
    - E_tot (kWh): array containing total energy vs. collector size
    - SWP_avg (L/kWh): array containing average efficiency vs. collector size (L/kWh)
    - W_vs_t (L/h): array with dimensions len(A)*len(t) containing all L/h vs. time traces
    - tplot (h): array containing times for W_vs_t, for plotting convenience

    Optional Outputs:
    - V_vs_t (L/h): if len(V_feed) !=0, output array with dimensions len(A)*len(t) containing all feed flow rates vs. 
        time traces
    - ratio_vs_t (1): if len(V_ratio) !=0, output array with dimensions len(A)*len(t) containing all flow rate ratios vs.
        time traces
    """

    # Give warning if dt not multiple of 5min
    if dt % 5 != 0:
        print(f'\nWarning: dt = {dt:0.1f} min not a multiple of 5 min, rounding dt to {int(np.around(dt/5, decimals=0)*5):0.0f} min')
    
    # Number of data points per time window
    npts = int(np.round(dt/5))

    # Solar collector sizes to test
    if A == []:
        A = np.arange(0.05,1+0.01,0.01) # m^2

    # Interpolated SWP vs. I_in for matching to each I_in vs. t pair
    fSWP = interp1d(P_in,SWP,kind=interp_kind) 

    # Interpolated flow rates and flow ratios if requested
    if len(V_feed) !=0:
        fVfeed = interp1d(P_in,V_feed,kind=interp_kind) 
    if len(V_ratio) !=0:
        fratio = interp1d(P_in,V_ratio,kind=interp_kind) 

    # Preallocate total production and average intensity array
    W = np.zeros(np.shape(A)) # L

    # Preallocate time traces arrays
    W_vs_t = np.zeros((len(A),len(t[0::npts])))
    V_vs_t = np.zeros((len(A),len(t[0::npts])))
    ratio_vs_t = np.zeros((len(A),len(t[0::npts])))
    
    # Loop over solar collector sizes
    for m in range(len(A)):

        # Get production over each time window
        W_dt = np.zeros(np.shape(t[0::npts])) # preallocate vector of production over dt
        V_dt = np.zeros(np.shape(t[0::npts])) # preallocate vector of flow rates over dt
        ratio_dt = np.zeros(np.shape(t[0::npts])) # preallocate vector of flow ratios over dt


        if dt == 5:
            for n in range(len(t)):
                if A[m]*I_in[n] >= np.amin(P_in) and A[m]*I_in[n] <= np.amax(P_in):
                    
                    # Store average SWP over dt
                    W_dt[n] = fSWP(A[m]*I_in[n]) * A[m]*I_in[n]/1000 * dt/60 # SWP*P_in*dt
                    
                    # Store flow rate and ratio if requested
                    if len(V_feed) != 0:
                        V_dt[n] = fVfeed(A[m]*I_in[n]) # V_feed
                    if len(V_ratio) != 0:
                        ratio_dt[n] = fratio(A[m]*I_in[n]) # V_ratio

        elif A[m]*I_in[n] > np.amax(P_in):
                    
                    # Store average SWP over dt
                    W_dt[n] = np.amax(SWP) * 1 * dt/60 # SWP*P_in*dt
                    
                    # Store flow rate and ratio if requested
                    if len(V_feed) != 0:
                        V_dt[n] = fVfeed(999*A[m]) # V_feed
                    if len(V_ratio) != 0:
                        ratio_dt[n] = fratio(999*A[m]) # V_ratio

        else:
            w = 0 # initialize window start index
            for n in range(len(t[0::npts])):

                if np.all(A[m]*I_in[w:w+npts] >= np.amin(P_in)):
                    W_dt[n] = np.trapz(fSWP(A[m]*I_in[w:w+npts]) * A[m]*I_in[w:w+npts]/1000 * dt/60, t[w:w+npts]) # SWP*P_in*dt
                
                    # Store flow rate and ratio if requested
                    if len(V_feed) != 0:
                        V_dt[n] = np.mean(fVfeed(A[m]*I_in[w:w+npts])) # V_feed
                    if len(V_ratio) != 0:
                        ratio_dt[n] = np.mean(fratio(A[m]*I_in[w:w+npts])) # V_ratio
                
                # Update window start index
                w = w + npts
 
        # Store time trace for current collector size
        W_vs_t[m,:] = W_dt/(dt/60)

        # Store flow rate and ratio time traces if requested
        V_vs_t[m,:] = V_dt
        ratio_vs_t[m,:] = ratio_dt
                
        # Record total production over time for current collector size
        W[m] = np.sum(W_dt)
    
    tplot = t[0::npts] # store discretized time array for plotting

    # Get average energy and efficiency vs collector size
    E_tot = np.zeros(np.shape(A))
    SWP_avg = np.zeros(np.shape(A))
    for m in range(len(A)):
        E_tot[m] = np.trapz(A[m]*I_in/1000,t) # kWh
        SWP_avg[m] = W[m]/E_tot[m]

    # Print results 
    if print_results == 1:
        print(f'Average solar intensity:  {np.mean(I_in):0.2f} W/m^2')
        print(f'Total energy:             {np.trapz(I_in/1000,t):0.2f} kWh/m^2/week')
        print(f'Average energy per day:   {np.trapz(I_in/1000,t)/7:0.2f} kWh/m^2/day')
        print(f'Time discretization:      {dt:0.0f} min')
        print(f'\n')
        print(f'Most efficient size:      {A[np.argmax(SWP_avg)]:0.2f} m^2')
        print(f'    Average SWP:          {np.amax(SWP_avg):0.2f} L/kWh/week')
        print(f'    Total production:     {W[np.argmax(SWP_avg)]:0.2f} L')
        print(f'    Total Energy:         {E_tot[np.argmax(SWP_avg)]:0.2f} kWh/week')
        print(f'\n')
        print(f'Most productive size:     {A[np.argmax(W)]:0.2f} m^2')
        print(f'    Average SWP:          {SWP_avg[np.argmax(W)]:0.2f} L/kWh/week')
        print(f'    Total production:     {np.amax(W):0.2f} L')
        print(f'    Total Energy:         {E_tot[np.argmax(W)]:0.2f} kWh/week')
        print(f'\n')

    if len(V_feed) !=0 and len(V_ratio) !=0: 
        return A, W, E_tot, SWP_avg, W_vs_t, tplot, V_vs_t, ratio_vs_t 
    elif len(V_ratio) !=0:
        return A, W, E_tot, SWP_avg, W_vs_t, tplot, ratio_vs_t
    elif len(V_feed) !=0:
        return A, W, E_tot, SWP_avg, W_vs_t, tplot, V_vs_t
    else:
        return A, W, E_tot, SWP_avg, W_vs_t, tplot
    









# Function to give performance for fixed flow rate and ratio for single collector size
def fixed_operation_performance(A,SWP,P_in,I_in,t,dt=5,print_results=1,V_feed=[],V_ratio=[],interp_kind='cubic'):

    """Inputs:
    - A (m^2): scalar value of single collector size
    - SWP (L/kWh): array containing maximum SWP below boiling of a given system vs. powers in P_in
    - P_in (W): powers at which the optimal SWP's were evaluated
    - I_in (W/m^2): solar intensities vs. times in t
    - t (h): times at which the intensities were measured

    Optional Inputs:
    - dt (min): time discretization/interval over which to evaluate average production rate (DEFAULT 5 min)
    - print_results: print aggregate results as text if true (DEFAULT True)
    - V_feed (L/h): array containing optimized feed flow rates vs. powers in P_in (DEFAULT EMPTY)
    - V_ratio (1): array containing optimized flow rate ratios vs. powers in P_in (DEFAULT EMPTY)"""

    # Give warning if dt not multiple of 5min
    if dt % 5 != 0:
        print(f'\nWarning: dt = {dt:0.1f} min not a multiple of 5 min, rounding dt to {int(np.around(dt/5, decimals=0)*5):0.0f} min')
    
    # Number of data points per time window
    npts = int(np.round(dt/5))

    # Interpolated SWP vs. I_in for matching to each I_in vs. t pair
    fSWP = interp1d(P_in,SWP,kind=interp_kind) 

    # Interpolated flow rates and flow ratios if requested
    if len(V_feed) !=0:
        fVfeed = interp1d(P_in,V_feed,kind=interp_kind) 
    if len(V_ratio) !=0:
        fratio = interp1d(P_in,V_ratio,kind=interp_kind) 

    # Preallocate time traces arrays
    W_vs_t = np.zeros(len(t[0::npts]))
    V_vs_t = np.zeros(len(t[0::npts]))
    ratio_vs_t = np.zeros(len(t[0::npts]))


    # Get production over each time window
    W_dt = np.zeros(np.shape(t[0::npts])) # preallocate vector of production over dt
    V_dt = np.zeros(np.shape(t[0::npts])) # preallocate vector of flow rates over dt
    ratio_dt = np.zeros(np.shape(t[0::npts])) # preallocate vector of flow ratios over dt

    if dt == 5:
        for n in range(len(t)):
            if A*I_in[n] >= np.amin(P_in):
                
                # Store average SWP over dt
                W_dt[n] = fSWP(A*I_in[n]) * A*I_in[n]/1000 * dt/60 # SWP*P_in*dt     

                # Store flow rate and ratio if requested
                if len(V_feed) != 0:
                    V_dt[n] = V_feed[0] 
                if len(V_ratio) != 0:
                    ratio_dt[n] = V_ratio[0]  

    else:
        w = 0 # initialize window start index
        for n in range(len(t[0::npts])):

            if np.all(A*I_in[w:w+npts] >= np.amin(P_in)):
                W_dt[n] = np.trapz(fSWP(A*I_in[w:w+npts]) * A*I_in[w:w+npts]/1000 * dt/60, t[w:w+npts]) # SWP*P_in*dt
                
                # Update window start index
                w = w + npts

                    # Store flow rate and ratio if requested
                if len(V_feed) != 0:
                    V_dt[n] = V_feed[0] 
                if len(V_ratio) != 0:
                    ratio_dt[n] = V_ratio[0]  
    

            
    # Store time trace for current collector size
    W_vs_t = W_dt/(dt/60)

    # Store flow rate and ratio time traces if requested
    V_vs_t = V_dt
    ratio_vs_t = ratio_dt
            
    # Record total production over time for current collector size
    W = np.sum(W_dt)
    
    tplot = t[0::npts] # store discretized time array for plotting

    # Get total energy and efficiency vs collector size
    E_tot = np.trapz(A*I_in/1000,t) # kWh
    SWP_avg = W/E_tot

    # Print results 
    if print_results == 1:
        print(f'Average solar intensity:  {np.mean(I_in):0.2f} W/m^2')
        print(f'Total energy:             {np.trapz(I_in/1000,t):0.2f} kWh/m^2/week')
        print(f'Average energy per day:   {np.trapz(I_in/1000,t)/7:0.2f} kWh/m^2/day')
        print(f'Time discretization:      {dt:0.0f} min')
        print(f'\n')
        print(f'Solar collector size:     {A:0.2f} m^2')
        print(f'    Average SWP:          {SWP_avg:0.2f} L/kWh/week')
        print(f'    Total production:     {W:0.2f} L')
        print(f'    Total Energy:         {E_tot:0.2f} kWh/week')
        print(f'\n')


    if len(V_feed) !=0 and len(V_ratio) !=0: 
        return W, E_tot, SWP_avg, W_vs_t, tplot, V_vs_t, ratio_vs_t 

    elif len(V_ratio) !=0:
        return W, E_tot, SWP_avg, W_vs_t, tplot, ratio_vs_t

    elif len(V_feed) !=0:
        return W, E_tot, SWP_avg, W_vs_t, tplot, V_vs_t

    else:
        return W, E_tot, SWP_avg, W_vs_t, tplot

