import spiceypy as spice
import numpy as np
from scipy.optimize import fsolve
import os
import time
import statistics
from math import sqrt
import matplotlib.pyplot as plt
from numpy.linalg import multi_dot
from typing import Optional, Tuple, Union
from quartic_solver import quartic_solver
import PICARD as picard

# Function transformation from et to tai and vice-versa
et2tai = np.vectorize(lambda t: spice.unitim(t, 'ET', 'TAI'), otypes=[float])
tai2et = np.vectorize(lambda t: spice.unitim(t, 'TAI', 'ET'), otypes=[float])
TWOPI = 2 * np.pi


def find_TCA_ET(mean_state_1: np.ndarray, mean_state_2: np.ndarray, TCA0_et: float) -> tuple:
    """
    Compute the Time of Closest Approach from the orthogonality between the 
    relative position and the relative velocity.
    
     INPUTS:
     - mean_state_1: primary object mean state in ECI reference frame at
                     TCA0_et
     - mean_state_2: secondary object mean state in ECI reference frame at
                     TCA0_et
      - TCA0_et     : TCA first guess (in ET)
 
      OUTPUT:
      - TCA_et_out: computed TCA (in ET)

      Author:  Andrea De Vittori, Politecnico di Milano, 21 March 2022
              e-mail: andrea.devittori@polimi.it
   """
    # Orbital elements
    elems_1 = spice.oscelt(mean_state_1, TCA0_et, 398600.4415)
    elems_2 = spice.oscelt(mean_state_2, TCA0_et, 398600.4415)

    # TCA0 in UTC
    TCA0_utc = TCA0_et - spice.deltet(TCA0_et, 'ET')

    # Set values
    prop_time_s_0 = -10
    iterations = 1
    max_iterations = 1000

    # Search for zeros untill convergence
    dot_prod_fun = lambda x: fsolve_fun(np.array(elems_1), np.array(elems_2), TCA0_utc, x) 
    flag_in = [1, 5]
    flag_out = 0
    # while not(flag_out in flag_in) & iterations<max_iterations:
    prop_time_s_0 = prop_time_s_0 + 10
    prop_time_s_out = fsolve(dot_prod_fun, prop_time_s_0, full_output=1, xtol=1e-5)
    flag_out = prop_time_s_out[2]
    iterations = iterations + 1
 
    # Raise an Error if fsolve failed
    # if  not(flag_out in flag_in):
    #     raise Exception('PoliMiSST:tcaIdentificaitonFailure:tcaComputationDidNotConverge TCA searching did not converge.\
    #      Consider modify the input TCA first guess.')

    # Time instants vect
    TCA_utc = TCA0_utc + prop_time_s_out[0]
    prop_time = TCA_utc - TCA0_utc
    TCA_et_out = TCA_utc + spice.deltet(TCA_utc,'UTC')
    mean_state_1_TCA = elems2state(elems_1, TCA_et_out, prop_time)
    mean_state_2_TCA = elems2state(elems_2, TCA_et_out, prop_time)
    return  TCA_et_out, TCA_utc, mean_state_1_TCA, mean_state_2_TCA

def fsolve_fun(elems_1: np.ndarray, elems_2: np.ndarray, TCA0_utc: float,
 prop_time_s: float) -> float:
    """
    Cost function to be minimized for the TCA estimation
    
    INPUTS:
    - elems_1     : primary osculating conic orbital elements at TCA0_et
    - elems_2     : secondary osculating conic orbital elements at TCA0_et
    - TCA0_utc    : reference TCA for elems_1 and elems_2
    - prop_time_s : forward or backward time window to minimize dr*dv

    OUTPUT:
    - drdv: dot product between dr and dv

    Author:  Andrea De Vittori, Politecnico di Milano, 21 March 2022
            e-mail: andrea.devittori@polimi.it
    """
    # % Propagation epoch in ET
    t_f_s = TCA0_utc + prop_time_s
    t_f_et = t_f_s + spice.deltet(t_f_s,'UTC')

    # % From orbital elements to states - primary
    x_1 = elems2state(elems_1, t_f_et, prop_time_s)

    # % From orbital elements to states - secondary
    x_2 = elems2state(elems_2, t_f_et, prop_time_s)

    # Compute relative position and velocity
    pos_rel = x_2[0:3]-x_1[0:3]
    vel_rel = x_2[3:6]-x_1[3:6]

    # Dot product between relative position and velocity
    drdv = np.dot(pos_rel, vel_rel)
    return drdv

def elems2state(elems: np.ndarray, t_f_et: float, prop_time: float) -> np.ndarray:
    """
     From conic orbital elements to state
     INPUTS:
     - elems       : osculating conic orbital elements at t_0_et
     - t_f_et      : final desired ephemeral time
     - prop_time   : forward or backward time in UTC span to link t_f_et to t_0_et
      OUTPUT:
      - state: state at t_f_et
      Author:  Andrea De Vittori, Politecnico di Milano, 21 March 2022
               e-mail: andrea.devittori@polimi.it
    """

    # Compute the mean anomaly to ridefine the state at  t_f_et
    a = elems[0]/(1-elems[1])
    n = np.sqrt(elems[7]/a**3)
    M = elems[5] + n*(prop_time)
    elems[5] = M
    elems[6] = t_f_et
    state = spice.conics(elems, t_f_et)
    return state

if __name__ == "__main__":

    # Lists definition
    err = []
    err_tot = []
    time_list = []
    name = []
    
    # Path for CDMs
    path = 'CDM_Python_Perturbed'
    CDMs = os.listdir('CDM_Python_Perturbed')

    # Load spice kernel
    spice.furnsh('latest_leapseconds.tls')

    # Parameters definition
    t0_et = 6e8
    mu = 398600.4415
    t0_utc = t0_et - spice.deltet(t0_et, 'ET')
    
    flag_quartic = 1

    # Method for picardlindelof
    coeff_extractor=picard.picardlindelof()
   
    # Loop over the CDMs 
    for i in CDMs:

        # Extract array for the Primary and Secondary states
        arr = np.loadtxt(path + '/' + i, delimiter=",")
        dim = len(arr)

        # Find the initial delta_t
        array = i.split("_")[1]
        delta_t = float(array.split(".txt")[0])
        name.append(delta_t)

        # Define the TCA guess
        t0_utc_guess = t0_utc - delta_t
        t0_et_guess = t0_utc_guess + spice.deltet(t0_utc_guess,'UTC')
        start = time.time()

        # ridefine err at each loop on i
        err = []
        for j in range(dim):

            # extract the primary and secondary states
            primary = arr[j][0:6]
            secondary = arr[j][6:]
            if flag_quartic:

                # Find TCA_et_new/TCA_UTC_new and Primary/Secondary state with the Quartic formula approximation
                a0,a1,a2,a3,a4=coeff_extractor.extractcoefficients(primary, secondary)
                z1,z2,z3,z4=quartic_solver(a0/a4,a1/a4,a2/a4,a3/a4)
                solutions=np.array([z1,z2,z3,z4])
                real_solutions= solutions[solutions.imag ==0]
                index = np.where(abs(real_solutions) == abs(real_solutions).min())
                computed_tcas =real_solutions[index].real
                elems_1 = spice.oscelt(primary, t0_et_guess, 398600.4415)
                elems_2 = spice.oscelt(secondary, t0_et_guess, 398600.4415)
                t0_utc_guess_new = t0_utc_guess + computed_tcas
                t0_et_guess_new = t0_utc_guess_new + spice.deltet(t0_utc_guess_new,'UTC')
                primary = elems2state(elems_1, t0_et_guess_new, computed_tcas)
                secondary = elems2state(elems_2,  t0_et_guess_new,  computed_tcas)


            # keplerian TCA finding
            TCA_out  = t0_utc_guess_new #find_TCA_ET(primary, secondary,  t0_et_guess_new)

            # Append the TCA error
            err.append(abs(float(TCA_out- t0_utc)))

        # Append the overall error and computational time
        err_tot.append(err)
        time_list.append(time.time()- start)

    # Plot the computational time vs TCA_guess
    fig = plt.figure()
    ax = plt.subplot()

    ax.scatter(name, time_list,  s=1000, color='black')
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)

    ax.set_xlabel('TCA guess [s]', fontsize=30)
    ax.set_ylabel('Computational time for 2170 cases [s]', fontsize=30)
    ax.grid()
    ax.set_axisbelow(True)

    # Plot the mean error vs TCA_guess
    fig = plt.figure()
    ax = plt.subplot()

    for i in range(len(CDMs)):
    
        array = CDMs[i].split("_")[1]
        
        delta_t = float(array.split(".txt")[0])
        print(delta_t) 
        ax.scatter(delta_t, np.mean(err_tot[i]),  s=500, color='black')
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)

    ax.set_xlabel('TCA guess [s]', fontsize=30)
    ax.set_ylabel('TCA mean estimation error [s]', fontsize=30)
    ax.grid()
    ax.set_axisbelow(True)

    plt.show()
