from cmath import isnan
from turtle import end_fill
import StochasticTaylorModel as stm
import numpy as np
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.linalg import pinvh
from scipy.integrate import odeint
import spiceypy as spice
import csv
from TCA_estimate import*
from collections import namedtuple
from quartic_solver import quartic_solver
import PICARD as picard
from mpl_toolkits.mplot3d import Axes3D
import random
from split_gmm import split_add
from sklearn import mixture
from sklearn.neighbors import KernelDensity

# Method for picardlindelof
coeff_extractor=picard.picardlindelof()

def PoC_Chan(x_p: np.ndarray, x_s: np.ndarray, P_ECI: np.ndarray, R: float, order_PoC: int) -> float:
    """
     PoC using the GMM formulation for short-term encounters
     INPUTS:
      - x_p         : Primary state in ECI rf. frame
      - x_s         : Secondary state in ECI rf. frame
      - P_ECI       : Summed Primary and Secondary Covariances in ECI rf.frame
      - R           : Summed primary and Secondary radii
      - order_PoC   : Chan's order expansion
    OUTPUTS:
      - P_m         : PoC
      Author:  Andrea De Vittori, Politecnico di Milano, 21 March 2022
               e-mail: andrea.devittori@polimi.it
    """

    # Primary and secondary state in ECI
    v_p = x_p[3:]
    v_s = x_s[3:]

    # Primary wrt the Secondary
    d_r_ECI = x_p[:3] - x_s[:3]

    # Primary wrt secondary in Bplane rf. frame
    z_vect = np.cross(v_s, v_p)/np.linalg.norm(np.cross(v_s, v_p))
    y_vect = (v_p - v_s)/np.linalg.norm(v_p - v_s)
    x_vect = np.cross(z_vect, y_vect)
    R_ECI2BP = np.vstack((z_vect, y_vect, x_vect))
    R_ECI2BP_2D =   np.array([[R_ECI2BP[0, 0], R_ECI2BP[0, 1], R_ECI2BP[0, 2]],\
         [ R_ECI2BP[2, 0], R_ECI2BP[2, 1], R_ECI2BP[2, 2]]])
    b_e = np.dot(R_ECI2BP_2D, d_r_ECI)

    # Covariance projection onto the Bplane rf. frame
    P_BP3 =  np.dot(np.dot(R_ECI2BP, P_ECI[:3, :3]), R_ECI2BP.T)
    P_BP2 =   np.array([[P_BP3[0, 0], P_BP3[0, 2]], [ P_BP3[2, 0], P_BP3[2, 2]]])
   
    # Definition of the parameter u and of the SMD v for the PoC formula
    u = R**2/(np.sqrt(P_BP2[0, 0]*P_BP2[1, 1])*np.sqrt(1-P_BP2[0, 1]**2/(P_BP2[0, 0]*P_BP2[1, 1])))
    v =  np.dot(np.dot(b_e, np.linalg.inv(P_BP2)), b_e)

    # Chan's PoC series initialization
    s_k = np.zeros([2,])
    S_k = np.zeros([2,])
    s_k[0] = np.exp(-u/2)
    S_k[0] = s_k[0]
    t_m = np.exp(-v/2)
    q_m = 1 - S_k[0]
    p_m = t_m*q_m
    P_m = p_m

    # Chan's PoC algorithm
    for m in range(1, order_PoC+1):
        s_k[1] = u/(2*m)*s_k[0]  
        S_k[1] = S_k[0] + s_k[1]
        t_m = v/(2*m)*t_m
        q_m = 1 -S_k[1]
        s_k[0] = s_k[1]
        S_k[0] = S_k[1]
        p_m = t_m*q_m
        P_m = P_m + p_m

    return  P_m

def PoC_GMM_ST(GMM_p: GaussianMixture, GMM_s: GaussianMixture, R: float, order_PoC: int) -> float:
    """
     PoC using the GMM formulation for short-term encounters
     INPUTS:
      - GMM_p       : Primary GMM
      - GMM_s       : Secondary GMM
      - R           : Summed Primary and Secondary radii
      - order_PoC   : Chan order expansion
     OUTPUT:
      - PoC_GMM     : PoC for the GMM
      Author:  Andrea De Vittori, Politecnico di Milano, 21 March 2022
               e-mail: andrea.devittori@polimi.it
    """
    # PoC_GMM initialization
    PoC_GMM = 0

    # Looping over the GMEs to compute the overall PoC_GMM
    for i in range(np.size(GMM_p.weights_)):
            PoC_GMM += GMM_p.weights_[i]*GMM_s.weights_[i]*\
                PoC_Chan(GMM_p.means_[i], GMM_s.means_[i], GMM_p.covariances_[i] + GMM_s.covariances_[i], R, order_PoC)
    return PoC_GMM

def PoC_ST(GMM_p: GaussianMixture, GMM_s: GaussianMixture, R: float, order_PoC: int,
 TCA_et: float, mu: float, flag_Kep: int) -> float:
    """
     PoC using the GMM formulation for short-term encounters
     INPUTS:
      - GMM_p       : Primary GMM
      - GMM_s       : Secondary GMM
      - R           : Summed Primary and Secondary radii
      - order_PoC   : Chan order expansion
      - TCA_et      : Primary and Secondary GMM ET time
      - mu          : Gravitational constant
      - flag_Kep    : FLag to activate the Keplerian refinement for TCA estimation
     OUTPUT:
      - PoC_GMM     : PoC for the GMM
      Author:  Andrea De Vittori, Politecnico di Milano, 21 March 2022
               e-mail: andrea.devittori@polimi.it
    """

    # Find number of mixtures for the Primary and Secondary
    mixtures_p = len(GMM_p.means_)
    mixtures_s = len(GMM_s.means_)

    # Define the struct for the Primary and Secondary pairs at TCA_et_new
    GMM_TCA = namedtuple("GMM_TCA" , "means_ covariances_ weights_ TCA_et_ TCA_UTC_")
    GMM_p_TCA = GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[])
    GMM_s_TCA = GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[])

    # Times definition
    TCA_UTC = TCA_et - spice.deltet(TCA_et,'ET')
    TCA_et_new = TCA_et

    # Looping over the Primary and Secondary GMMs
    for i in range(mixtures_p):
        for j in range(mixtures_s):

            # Append weights if the input GMMs to GMM_p_TCA and GMM_s_TCA
            GMM_p_TCA.weights_.append(GMM_p.weights_[i])
            GMM_s_TCA.weights_.append(GMM_s.weights_[j])

            # Primary and secondary state definition in ECI rf. frame
            primary = GMM_p.means_[i]
            secondary =  GMM_s.means_[j]
            try:
            # Find TCA_et_new/TCA_UTC_new and Primary/Secondary state with the Quartic formula approximation
                a0,a1,a2,a3,a4=coeff_extractor.extractcoefficients(primary, secondary)
                z1,z2,z3,z4=quartic_solver(a0/a4,a1/a4,a2/a4,a3/a4)
                solutions=np.array([z1,z2,z3,z4])
                real_solutions= solutions[solutions.imag ==0]
                if real_solutions.any():
                    index = np.where(abs(real_solutions) == abs(real_solutions).min())
                    computed_tcas = real_solutions[index].real
                    elems_1 = spice.oscelt(GMM_p.means_[i], TCA_et, 398600.4415)
                    elems_2 = spice.oscelt(GMM_s.means_[j], TCA_et, 398600.4415)
                    TCA_UTC_new = TCA_UTC + computed_tcas
                    TCA_et_new = TCA_UTC_new + spice.deltet(TCA_UTC_new,'UTC')
                    primary = elems2state(elems_1, TCA_et_new,  computed_tcas)
                    secondary = elems2state(elems_2, TCA_et_new,   computed_tcas)
                if (flag_Kep and real_solutions.any()) or (not real_solutions.any()):
                    TCA_out = find_TCA_ET(primary, secondary, TCA_et_new)
                    TCA_et_new = TCA_out[0]
                    TCA_UTC_new = TCA_out[1]
                    primary = TCA_out[2]
                    secondary = TCA_out[3]
            except:
                TCA_out = find_TCA_ET(primary, secondary, TCA_et_new)
                TCA_et_new = TCA_out[0]
                TCA_UTC_new = TCA_out[1]
                primary = TCA_out[2]
                secondary = TCA_out[3]
            # Refinement of the Quartic formula approxiamtion adopting a keplerian TCA finding


            # Append to GMM_p_TCA and GMM_s_TCA the state at TCA and the correspondig TCA_et_new/TCA_UTC_new
            GMM_p_TCA.means_.append(primary)
            GMM_s_TCA.means_.append(secondary)
            GMM_p_TCA.TCA_et_.append(TCA_et_new)
            GMM_p_TCA.TCA_et_.append(TCA_UTC_new)
            GMM_s_TCA.TCA_et_.append(TCA_et_new)
            GMM_s_TCA.TCA_et_.append(TCA_UTC_new)

            # Append the propagated covariance using an analytic STM
            GMM_p_TCA.covariances_.append(propagate_covariance(GMM_p.means_[i], 
                GMM_p.covariances_[i], TCA_et,  TCA_et_new, mu))
            GMM_s_TCA.covariances_.append(propagate_covariance(GMM_s.means_[j], 
                GMM_s.covariances_[j], TCA_et,  TCA_et_new, mu))

    # PoC computation for short-term encounters
    PoC_out = PoC_GMM_ST(GMM_p_TCA, GMM_s_TCA, R, order_PoC)
    return PoC_out

if __name__ == "__main__":

    # Loading latest_leapsecond kernels
    spice.furnsh('latest_leapseconds.tls')

    # Initialization
    flag = np.array([3, 1, 0, 0, 0, 0, 1])
    flag_refinement = 0
    n_splits_p = 9 # Maximum number of splits for the primary
    n_splits_s = 9 # Maximum number of splits for the secondary
    n_MC = 10000 # Montecarlo sampling number
    PoC_order = 6
    t_UTC_0 = 6e8
    mu = 398600.4415
    np.random.seed(1)
    GMM_TCA = namedtuple("GMM_TCA" , "means_ covariances_ weights_ TCA_et_ TCA_UTC_")
    
    if flag[0] == 0:
        prop_time = 280800

        R_sum = 15e-3

        PoC_true = 0.100846420
        
        x_p_0= np.array([-33552760.2666400000, -23727886.0120850000, 0.0000000000,
       -1828.9654659653, 2534.1298970058, 0.0000000000])/1000# [km/s]

        x_s_0 = np.array([-33161707.8066380000, -24056602.5363940000, 122462.5305735300,
           -1859.4197012053, 2523.7257583343, 6.9196063841])/1000

        x_p_TCA= np.array([153951.4752631000, 41874153.9947520000, 0.0000000000,
           3066.8746235984, -11.4110245828, 0.00000000009])/1000# [km/s]

        x_s_TCA = np.array([153951.9732739700, 41874156.7446250000, 2.7520743270,
                3066.8646233948, -0.0449985917, -11.3560272115])/1000

        P_p_0 = np.array([[0.0571247004665740000000, -0.0237271923593760000000, 0.0000000000000000000000,   0.0000000000000000000000,   0.0000000000000000000000,   0.0000000000000000000000], 
                [-0.0237271923593760000000, 0.0728752995334260000000, 0.0000000000000000000000,   0.0000000000000000000000,   0.0000000000000000000000,   0.0000000000000000000000],
                [0.0000000000000000000000, 0.0000000000000000000000, 0.0400000000000000000000,   0.0000000000000000000000,   0.0000000000000000000000,   0.0000000000000000000000], 
                [0.0000000000000000000000, 0.0000000000000000000000, 0.0000000000000000000000,   0.0000000100000000000000,   0.0000000000000000000000,   0.0000000000000000000000],    
                [0.0000000000000000000000, 0.0000000000000000000000, 0.0000000000000000000000,   0.0000000000000000000000 ,  0.0000000100000000000000,   0.0000000000000000000000],    
                [0.000000000000000000000, 0.0000000000000000000000, 0.0000000000000000000000  , 0.0000000000000000000000,   0.0000000000000000000000,   0.0000000100000000000000]])/1000**2

        P_s_0 = np.array([[0.0575921124921660000000, -0.0238771630800800000000 ,  -0.0000654669270369290000   ,0.0000000000000000000000,   0.0000000000000000000000 ,  0.0000000000000000000000], 
                            [-0.0238771630800800000000, 0.0724076438805530000000, 0.0000888559855394610000  , 0.0000000000000000000000 ,  0.0000000000000000000000  , 0.0000000000000000000000],
                            [-0.0000654669270369290000, 0.0000888559855394610000, 0.0400002436272810000000 ,  0.0000000000000000000000 ,  0.0000000000000000000000 ,  0.0000000000000000000000], 
                            [0.0000000000000000000000,0.0000000000000000000000, 0.0000000000000000000000 ,  0.0000000100000000000000  , 0.0000000000000000000000  , 0.0000000000000000000000],    
                            [0.0000000000000000000000, 0.0000000000000000000000,0.0000000000000000000000 ,  0.0000000000000000000000 ,  0.0000000100000000000000  , 0.0000000000000000000000],    
                            [0.0000000000000000000000,0.0000000000000000000000,0.0000000000000000000000,   0.0000000000000000000000,   0.0000000000000000000000,   0.0000000100000000000000]])/1000**2

        P_p_TCA = np.array([[6494.0687699936000000000000, -376.2174293862200000000000, 0.0000000000000000000000,   0.0159835090297390000000,  -0.4942613967814000000000,   0.0000000000000000000000], 
                    [-376.2174293862200000000000, 22.5686372188410000000000, 0.0000000000000000000000,  -0.0009881643469155700000 ,  0.0285754945542590000000 ,  0.0000000000000000000000],
                    [0.0000000000000000000000, 0.0000000000000000000000, 1.2050405732107000000000 ,  0.0000000000000000000000 ,  0.0000000000000000000000,  -0.0000607087721201980000], 
                    [0.0159835090297390000000, -0.0009881643469155700000, 0.0000000000000000000000 ,  0.0000000443435166923780,  -0.0000012117808283202000 ,  0.0000000000000000000000],    
                    [-0.4942613967814000000000, 0.0285754945542590000000, 0.0000000000000000000000 , -0.0000012117808283202000 ,  0.0000376229959779560000  , 0.0000000000000000000000],    
                    [0.0000000000000000000000, 0.0000000000000000000000, -0.0000607087721201980000 ,  0.0000000000000000000000 ,  0.0000000000000000000000,   0.0000000033903879281480]])/1000**2


        P_s_TCA = np.array([[6539.7159476886000000000000, -354.4827622708000000000000 ,  -24.2158339268500000000000,   0.0161364407930240000000 , -0.4975675313434400000000,  -0.0000667805505728000000], 
                    [-354.4827622708000000000000, 19.9876841613490000000000, 1.3128407317001000000000,  -0.0009366299987872100000   ,0.0269117100210450000000 ,  0.0000038371255372188000],
                    [-24.2158339268500000000000,1.3128407317001000000000, 1.2674955658292000000000 , -0.0000599870075272560000 ,  0.0018427685478909000000,  -0.0000602391913489130000], 
                    [0.0161364407930240000000, -0.0009366299987872100000, -0.0000599870075272560000,   0.0000000447836341841990,  -0.0000012229871227354000 , -0.0000000001697954617528],    
                    [-0.4975675313434400000000, 0.0269117100210450000000, 0.0018427685478909000000,  -0.0000012229871227354000  , 0.0000378619099092410000 ,  0.0000000050463685349191],    
                    [-0.0000667805505728000000, 0.0000038371255372188000, -0.0000602391913489130000  ,-0.0000000001697954617528  , 0.0000000050463685349191  , 0.0000000034465981691675]])/1000**2

    elif flag[0] == 1:

        prop_time = 172800

        R_sum = 10e-3

        PoC_true = 0.004300500

        x_p_0= np.array([-6345736.9319327000,-1876218.4567378000,-1876218.4567378000,
                  2936.7099630585,-4966.2630709530,-4966.2630709530])/1000# [km/s]

        x_s_0 = np.array([-6345839.9792600000,-1876012.2715538000,-1875945.3383511000,
                  2936.4524265637,-4966.4623670708,-4966.2765632153])/1000
        
        x_p_TCA= np.array([6877715.6342470000, 53834.2140394360, 53834.2140394360,
                 -84.2628278054, 5382.5970952289, 5382.5970952289])/1000# [km/s]

        x_s_TCA = np.array([6877716.6344015000, 53835.2141405120, 53836.2137729820,
                 -84.1628151363, 5382.6970892007,5382.4970647006])/1000

        P_p_0 = np.array([[4.7440894789163000000000, -1.2583279067770000000000,   -1.2583279067770000000000,   0.0000000000000000000000,   0.0000000000000000000000,   0.0000000000000000000000], 
                [-1.2583279067770000000000,6.1279552605419000000000,2.1279552605419000000000 ,  0.0000000000000000000000,   0.0000000000000000000000 ,  0.0000000000000000000000],
                [-1.2583279067770000000000,2.1279552605419000000000,6.1279552605419000000000,   0.0000000000000000000000,   0.0000000000000000000000,   0.0000000000000000000000], 
                [0.0000000000000000000000,0.0000000000000000000000,0.0000000000000000000000,   0.0000010000000000000000,   0.0000000000000000000000,   0.0000000000000000000000],    
                [0.0000000000000000000000,0.0000000000000000000000,0.0000000000000000000000,   0.0000000000000000000001,   0.0000010000000000000000,  -0.0000000000000000000001],    
                [0.0000000000000000000000,0.0000000000000000000000,0.0000000000000000000000,   0.0000000000000000000001,  -0.0000000000000000000001 ,  0.0000010000000000000000]])/1000**2

        P_s_0 = np.array([[4.7439512624715000000000,  -1.2582550000046000000000, -1.2582079265320000000000,  0.0000000000000000000000,   0.0000000000000000000000, 0.0000000000000000000000], 
                            [-1.2582550000046000000000, 6.1281039832865000000000,  2.1280243672750000000000,   0.0000000000000000000000,   0.0000000000000000000000, 0.0000000000000000000000],
                            [-1.2582079265320000000000, 2.1280243672750000000000,  6.1279447542420000000000,   0.0000000000000000000000,   0.0000000000000000000000, 0.0000000000000000000000], 
                            [0.0000000000000000000000,  0.0000000000000000000000,  0.0000000000000000000000,   0.0000010000000000000000,   0.0000000000000000000000, 0.0000000000000000000000],    
                            [0.0000000000000000000000,  0.0000000000000000000000,  0.0000000000000000000000,   0.0000000000000000000000,   0.0000010000000000000000, -0.0000000000000000000002],    
                            [0.0000000000000000000000,  0.0000000000000000000000,  0.0000000000000000000000,   0.0000000000000000000000,  -0.0000000000000000000002, 0.0000010000000000000000]])/1000**2

        P_p_TCA = np.array([[430.02827939777, -18393.891381385,   -18393.891381385, 28.811775614851,     0.14985494859082,   0.14985494859082], 
                    [-18393.891381385, 790190.51867865, 790186.9676547,  -1237.8765173077,    -6.5012185382348,   -6.5024449244912],
                    [-18393.891381385, 790186.9676547,  790190.51867865, -1237.8765173077,   -6.5024449244912,   -6.5012185382348], 
                    [28.811775614851, -1237.8765173077,  -1237.8765173077,1.9392157223809,   0.010188218207003,   0.010188218207003],    
                    [0.14985494859082, -6.5012185382348,  -6.5024449244912, 0.010188218207003,    0.000055474681548641,   0.000053924699534112],    
                    [0.14985494859082, -6.5024449244912, -6.5012185382348, 0.010188218207003,  0.000053924699534112,    0.000055477431358815]])/1000**2


        P_s_TCA = np.array([[429.35284572627, -18379.895353487,   -18379.212284332, 28.78930421662,     0.14973961436546,   0.14974660730687], 
                    [-18379.895353487, 790234.93895821,  790202.01968289, -1237.9226383171,   -6.5015746157359,   -6.503099317208],
                    [-18379.212284332, 790202.01968289,   790176.20377729, -1237.8766321984,   -6.5025593138425,   -6.5016312695197], 
                    [28.78930421662, -1237.922638317, -1237.8766321984,  1.9392512139259,    0.010188582861453,   0.010189050148708],    
                    [0.14973961436546, -6.5015746157359,  -6.5025593138425, 0.010188582861453,    0.000055477431358815,    0.00005392995296381],    
                    [0.14974660730687, -6.503099317208,  -6.5016312695197, 0.010189050148708,  0.00005392995296381,    0.000055477431358815]])/1000**2
    elif flag[0] == 2:

        R_sum = 2.97e-2

        PoC_true =  1.383503389181375e-01

        x_p_TCA= np.array([    2.330521851751370e+00, -1.103704510502010e+03 , 7.105887642997180e+03,
    -7.442862828717730e+00, -6.137347436526600e-04, 3.951361392933490e-03])

        x_s_TCA = np.array([    2.333465506263320e+00,-1.103671212478360e+03,7.105914958099040e+03,
     7.353740487126320e+00, -1.142814049765360e+00, -1.982472259113770e-01])

        P_p_TCA = np.array([[  1.777981361085467e-02,    -1.315989244494694e-04 ,    2.392444609963430e-04, 0,0,0],
                            [-1.315989244494694e-04    , 2.806338650876504e-05  ,  -3.365999545058470e-05, 0,0,0],
                            [2.392444609963430e-04 ,   -3.365999545058469e-05 ,    8.411224349521375e-05, 0,0,0],
                            [0,0,0,0,0,0],
                            [0,0,0,0,0,0],
                            [0,0,0,0,0,0]])


        P_s_TCA = np.array([[     8.004472543672787e-01  ,  -1.232355787216603e-01    ,-2.135367751841215e-02, 0,0,0],
                            [-1.232355787216603e-01   ,  1.921491693112587e-02   ,  3.301796191822024e-03,0,0,0],
                            [-2.135367751841215e-02   ,  3.301796191822024e-03  ,   1.213456190606059e-03,0,0,0],
                            [0,0,0,0,0,0],
                            [0,0,0,0,0,0],
                            [0,0,0,0,0,0]])
    elif flag[0] == 3:
        R_sum = 0.005000000000000*2 
        PoC_true = 6.090669795698510e-04
        x_p_TCA= np.array([  4200.71593969004,
-4979.91690172379,
-1668.73897456795,
2.68467953119580,
4.21791276477550,
-5.85085409673824])

        x_s_TCA = np.array([    4200.96297696341,
-4979.83531885185,
-1668.55962545501,
-2.88802435960042,
-4.33899807328328,
5.71743802503703])

        P_p_TCA = np.array([[0.00659699302384435,	0.0104094522144037,	-0.0143868158958461,	-1.35576239547068e-05,	1.60571647366398e-05,	5.35412021511665e-06],
[0.0104094522144037,	0.0164581864535499,	-0.0227573607787527,	-2.14315213846218e-05,	2.53952516441365e-05,	8.45696521903618e-06],
[-0.0143868158958461,	-0.0227573607787527,	0.0315152007877405,	2.96604305211574e-05,	-3.51396289797674e-05,	-1.17110490781476e-05],
[-1.35576239547068e-05,	-2.14315213846218e-05,	2.96604305211574e-05,	2.79483350131319e-08,	-3.30747765832376e-08,	-1.10102149522667e-08],
[1.60571647366398e-05,	2.53952516441365e-05,	-3.51396289797674e-05,	-3.30747765832376e-08,	3.92085111849983e-08,	1.30681551276307e-08],
[5.35412021511665e-06,	8.45696521903618e-06,	-1.17110490781476e-05,	-1.10102149522667e-08,	1.30681551276307e-08,	4.37169828135763e-09]])


        P_s_TCA = np.array([[0.562763133844756,	0.821368403185235,	-1.06744933246535,	0.00103989059858363,	-0.00122417165666914,	-0.000413195718782975],
[0.821368403185235,	1.22711388212953,	-1.60359532905332,	0.00154648287251050,	-0.00183533342027384,	-0.000610502685985866],
[-1.06744933246535,	-1.60359532905332,	2.14220317941386,	-0.00204477850856149,	0.00242174598218732,	0.000815531270990066],
[0.00103989059858363,	0.00154648287251050,	-0.00204477850856149,	1.98056996073545e-06,	-2.32087977576463e-06,	-7.73034137142535e-07],
[-0.00122417165666914,	-0.00183533342027384,	0.00242174598218732,	-2.32087977576463e-06,	2.76589485247476e-06,	9.29669467431308e-07],
[-0.000413195718782975,	-0.000610502685985866,	0.000815531270990066,	-7.73034137142535e-07,	9.29669467431308e-07,	3.25581833463626e-07]])

    elif flag[0] == 4:
        R_sum = 0.005000000000000*2 
        PoC_true = 0.005029280444607
        x_p_TCA= np.array([ -42896.7517778083,
-1875.59187201092,
-354.626075876183,
0.0591020703965643,
-1.70879100376313,
-0.176045674771793])

        x_s_TCA = np.array([-42896.7616990771,
-1875.58784660205,
-354.658625442134,
-0.150069757746341,
3.01985124101414,
0.472502720147411])

        P_p_TCA = np.array([[0.0276163693576352,	-0.615289243035881,	-0.0643588925890301,	7.89615644555697e-05,	4.55925093113454e-06,	5.34914212344894e-07],
[-0.615289243035881,	13.8635109963118,	1.42749144315257,	-0.00177901362093960,	-0.000103583860831128, -1.20158176960245e-05],
[-0.0643588925890301,	1.42749144315257,	0.176934057707096,	-0.000183111554754670,	-1.05716616014874e-05,	-1.14914484889313e-06],
[7.89615644555697e-05,	-0.00177901362093960,	-0.000183111554754670,	2.28535237521285e-07,	1.33075917540813e-08,	1.54837094784773e-09],
[4.55925093113454e-06,	-0.000103583860831128,	-1.05716616014874e-05,	1.33075917540813e-08,	7.81294684700552e-10,	7.52964400259622e-11],
[5.34914212344894e-07,	-1.20158176960245e-05,	-1.14914484889313e-06,	1.54837094784773e-09,	7.52964400259622e-11,	1.47171387501996e-10]])


        P_s_TCA = np.array([[0.000730181891745581,	-0.0109450910690001,	-0.00214346821546025,	-7.90396656695215e-07,	-2.32328184423784e-08,	5.20002606428932e-08],
[-0.0109450910690001,	0.404181028224225,	0.0562690363782541,	2.87729493074076e-05,	1.52969002889315e-06,	8.06479582657675e-07],
[-0.00214346821546025,	0.0562690363782541,	0.0651911904677833,	4.70690048457187e-06,	2.31137116754839e-07,	9.44524789368751e-08],
[-7.90396656695216e-07,	2.87729493074076e-05,	4.70690048457187e-06,	2.05904128071157e-09,	1.09808245148057e-10,	5.15894125172733e-11],
[-2.32328184423784e-08,	1.52969002889315e-06,	2.31137116754839e-07,	1.09808245148057e-10,	1.33567027887199e-11,	-3.68061601888716e-11],
[5.20002606428932e-08,	8.06479582657675e-07,	9.44524789368751e-08,	5.15894125172732e-11,	-3.68061601888717e-11,	2.88529665878764e-10]])



    if flag[1]:
        GMM_p =[]
        GMM_s =[]
        # Define the struct for the Primary and Secondary pairs at TCA_et_new
        GMM_p = GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[])
        GMM_s = GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[])
        GMM_p.means_.append(x_p_TCA)
        GMM_p.covariances_.append(P_p_TCA)
        GMM_p.weights_.append(1)
        GMM_s.means_.append(x_s_TCA)
        GMM_s.covariances_.append(P_s_TCA)
        GMM_s.weights_.append(1)
        # Test Chan at TCA
        PoC = PoC_GMM_ST(GMM_p, GMM_s, R_sum, PoC_order)
        plt.figure()
        ax = plt.subplot(111, projection='3d')
        dx = np.ones(1)
        dy = np.ones(1)
        xpos = np.array([0]).flatten()
        ypos = np.array([0]).flatten()
        zpos = np.zeros(1)
        dz = PoC.flatten()
        X = np.arange(-0.5, 2.5, 1)
        Y = np.arange(-0.5, 2.5, 1)
        X, Y = np.meshgrid(X, Y)
        Z = np.ones((3, 3))*PoC_true
        cmap = plt.cm.get_cmap('jet') # Get desired colormap - you can change this!
        max_height = np.max(dz)   # get range of colorbars so we can normalize
        min_height = np.min(dz)
            # scale each z to [0,1], and get their rgb values
        rgba = [cmap((k-min_height)/max_height) for k in dz] 
        surf = ax.plot_surface(X, Y, Z, color='green',
                    linewidth=0, antialiased=False, alpha = 0.5)
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz,  color=rgba, zsort='average')
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.tick_params(axis='z', labelsize=10)
        ax.set_xlabel('N Primary', fontsize=10)
        ax.set_ylabel('N secondary', fontsize=10)
        ax.set_zlabel('Probability', fontsize=10)

    # if flag[2:].any():
            # Test Chan with GMMa and MC simulation
            # x_p_MC_0 = np.random.multivariate_normal(x_p_0, P_p_0, n_MC).T
            # x_s_MC_0 = np.random.multivariate_normal(x_s_0, P_s_0, n_MC).T
            # t_ET_0 = t_UTC_0 + spice.deltet(t_UTC_0,'UTC')
            # t_UTC_TCA = t_UTC_0 + prop_time
           # t_et_TCA =  t_UTC_TCA + spice.deltet( t_UTC_TCA,'UTC')
    #         x_p_MC_TCA = np.zeros((6, n_MC))
    #         x_s_MC_TCA = np.zeros((6, n_MC))
    #         for i in range(n_MC):
    #             elems = spice.oscelt(x_p_MC_0[:, i], t_ET_0, mu)
    #             x_p_MC_TCA[:, i] = elems2state(elems, t_et_TCA, prop_time)
    #             elems = spice.oscelt(x_s_MC_0[:, i], t_ET_0, mu)
    #             x_s_MC_TCA[:, i] = elems2state(elems, t_et_TCA, prop_time)

    #         P_p_TCAv = np.cov(x_p_MC_TCA)
    #         P_s_TCAv = np.cov(x_s_MC_TCA)
    #         x_p_TCAv = np.mean(x_p_MC_TCA, axis = 1)
    #         x_s_TCAv = np.mean(x_s_MC_TCA, axis = 1)

    if flag[2]:
        GMM_p_TCA = []
        GMM_s_TCA = []
        PoC_out =np.zeros((n_splits_p, n_splits_s))
        #   GMM for the primary object
        for i in range(1,n_splits_p+1):
            # GMM_p_TCA.append(mixture.BayesianGaussianMixture(n_components=i).fit(x_p_MC_TCA.T))
            a = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(x_p_MC_TCA.T)
            # GMM_p_TCA.append(GaussianMixture(n_components=i, tol = 1e-17).fit(x_p_MC_TCA.T))
        for i in range(1,n_splits_s+1):
            # GMM for the primary object
            GMM_s_TCA.append(mixture.BayesianGaussianMixture(n_components=i).fit(x_s_MC_TCA.T))
            # GMM_s_TCA.append(GaussianMixture(n_components=i, tol = 1e-17).fit(x_s_MC_TCA.T))
        for i in range(0,n_splits_p):
            for j in range(0,n_splits_s):
                PoC_out[i, j] = PoC_ST(GMM_p_TCA[i],  GMM_s_TCA[j], R_sum,  PoC_order,  t_et_TCA, mu, flag_refinement)

        plt.figure()
        ax = plt.subplot(111, projection='3d')
        dx = np.ones(n_splits_p*n_splits_s)
        dy = np.ones(n_splits_p*n_splits_s)
        xpos = np.array([np.arange(0, n_splits_p, 1)]*n_splits_s).flatten()
        ypos = np.array([[i]*n_splits_p for i in range(n_splits_s)]).flatten()
        zpos = np.zeros(n_splits_p*n_splits_s)
        dz = PoC_out.flatten()
        X = np.arange(0, n_splits_p+1, 1)
        Y = np.arange(0, n_splits_p+1, 1)
        X, Y = np.meshgrid(X, Y)
        Z = np.ones((n_splits_p+1,n_splits_p+1))*PoC_true
        cmap = plt.cm.get_cmap('jet') # Get desired colormap - you can change this!
        max_height = np.max(dz)   # get range of colorbars so we can normalize
        min_height = np.min(dz)
            # scale each z to [0,1], and get their rgb values
        rgba = [cmap((k-min_height)/max_height) for k in dz] 
        surf = ax.plot_surface(X, Y, Z, color='green',
                    linewidth=0, antialiased=False, alpha = 0.5)
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz,  color=rgba, zsort='average')
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.tick_params(axis='z', labelsize=10)
        ax.set_xlabel('N Primary', fontsize=10)
        ax.set_ylabel('N secondary', fontsize=10)
        ax.set_zlabel('Probability', fontsize=10)    
            
    if flag[3]:
    # Case Vittaldev + GMM

        # Define the struct for the Primary and Secondary pairs at TCA_et_new
        GMM_p_TCAv = []
        GMM_s_TCAv = []
        w_p,v_p=np.linalg.eig(P_p_TCAv)
        index_p = np.where(w_p == w_p.max())
        GMM_pv = []
        n=0
        for i in range(3,n_splits_p+1, 2):
            GMM_pv.append(GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[]))
            xi,Pi,wi = split_add(x_p_TCAv, P_p_TCAv, np.squeeze(v_p[:,index_p]), i)
            GMM_p_TCAv.append(GaussianMixture(n_components=i, tol = 1e-17, means_init=xi, weights_init=wi).fit(x_p_MC_TCA.T))
            for j in range(i):
                GMM_pv[n].means_.append(xi[j, :])
                GMM_pv[n].covariances_.append(Pi[j, :, :])
                GMM_pv[n].weights_.append(wi[j])
            n +=1

        w_s,v_s=np.linalg.eig(P_s_TCAv)
        index_s = np.where(w_s == w_s.max())
        GMM_sv =[]
        n = 0
        for i in range(3,n_splits_s+1, 2):
            GMM_sv.append(GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[]))
            xi,Pi,wi = split_add(x_s_TCAv, P_s_TCAv, np.squeeze(v_s[:,index_s]), i)
            GMM_s_TCAv.append(GaussianMixture(n_components=i, tol = 1e-17, means_init=xi, weights_init=wi).fit(x_s_MC_TCA.T))
            for j in range(i):
                GMM_sv[n].means_.append(xi[j, :])
                GMM_sv[n].covariances_.append(Pi[j, :, :])
                GMM_sv[n].weights_.append(wi[j])
            n +=1

        n_el_p = len(GMM_p_TCAv)
        n_el_s = len(GMM_s_TCAv)
        PoC_outv = np.zeros((n_el_p*n_el_s, ))
        n = 0
        for j in range(0,n_el_s):
            for i in range(0,n_el_p):
                PoC_outv[n] = PoC_ST(GMM_p_TCAv[i],  GMM_s_TCAv[j], R_sum,  PoC_order,  t_et_TCA, mu, flag_refinement)
                n += 1

        plt.figure()
        ax = plt.subplot(111, projection='3d')
        dx = np.ones(n_el_p*n_el_s)
        dy = np.ones(n_el_p*n_el_p)
        xpos = np.array([np.arange(2, n_splits_p+1, 2)]*n_el_s).flatten()
        ypos = np.array([[i]*n_el_p for i in range(2, n_splits_s+1, 2)]).flatten()
        zpos = np.zeros(n_el_p*n_el_s)
        dz = PoC_outv
        X = np.arange(2, 2+2*n_el_p)
        Y = np.arange(2, 2+2*n_el_s)
        X, Y = np.meshgrid(X, Y)
        Z = np.ones((2*n_el_s,2*n_el_p))*PoC_true
        cmap = plt.cm.get_cmap('jet') # Get desired colormap - you can change this!
        max_height = np.max(dz)   # get range of colorbars so we can normalize
        min_height = np.min(dz)
            # scale each z to [0,1], and get their rgb values
        rgba = [cmap((k-min_height)/max_height) for k in dz] 
        surf = ax.plot_surface(X, Y, Z, color='green',
                        linewidth=0, antialiased=False, alpha = 0.5)
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz,  color=rgba, zsort='average')
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.tick_params(axis='z', labelsize=10)
        ax.set_xlabel('N Primary', fontsize=10)
        ax.set_ylabel('N secondary', fontsize=10)
        ax.set_zlabel('Probability', fontsize=10)

    if flag[4]:
        GMM_pv = []
        w_p,v_p=np.linalg.eig(P_p_TCAv)
        index_p = np.where(w_p == w_p.max())
        n=0

        for i in range(3,n_splits_p+1, 2):
            GMM_pv.append(GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[]))
            xi,Pi,wi = split_add(x_p_TCAv, P_p_TCAv, np.squeeze(v_p[:,index_p]), i)
            for j in range(i):
                GMM_pv[n].means_.append(xi[j, :])
                GMM_pv[n].covariances_.append(Pi[j, :, :])
                GMM_pv[n].weights_.append(wi[j])
            n +=1

        w_s,v_s=np.linalg.eig(P_s_TCAv)
        index_s = np.where(w_s == w_s.max())
        GMM_sv =[]
        n = 0
        for i in range(3,n_splits_s+1, 2):
            GMM_sv.append(GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[]))
            xi,Pi,wi = split_add(x_s_TCAv, P_s_TCAv, np.squeeze(v_s[:,index_s]), i)
            for j in range(i):
                GMM_sv[n].means_.append(xi[j, :])
                GMM_sv[n].covariances_.append(Pi[j, :, :])
                GMM_sv[n].weights_.append(wi[j])
            n +=1

        n = 0
        n_el_p = len(GMM_pv)
        n_el_s = len(GMM_sv)
        PoC_outvv = np.zeros((n_el_p*n_el_s, ))
        for j in range(0,n_el_s):
            for i in range(0,n_el_p):
                PoC_outvv[n] = PoC_ST(GMM_pv[i],  GMM_sv[j], R_sum,  PoC_order,  t_et_TCA, mu, flag_refinement)
                n += 1

        plt.figure()
        ax = plt.subplot(111, projection='3d')
        dx = np.ones(n_el_p*n_el_s)
        dy = np.ones(n_el_p*n_el_p)
        xpos = np.array([np.arange(2, n_splits_p+1, 2)]*n_el_s).flatten()
        ypos = np.array([[i]*n_el_p for i in range(2, n_splits_s+1, 2)]).flatten()
        zpos = np.zeros(n_el_p*n_el_s)
        dz = PoC_outvv
        X = np.arange(2, 2+2*n_el_p)
        Y = np.arange(2, 2+2*n_el_s)
        X, Y = np.meshgrid(X, Y)
        Z = np.ones((2*n_el_s,2*n_el_p))*PoC_true
        cmap = plt.cm.get_cmap('jet') # Get desired colormap - you can change this!
        max_height = np.max(dz)   # get range of colorbars so we can normalize
        min_height = np.min(dz)
            # scale each z to [0,1], and get their rgb values
        rgba = [cmap((k-min_height)/max_height) for k in dz] 
        surf = ax.plot_surface(X, Y, Z, color='green',
                        linewidth=0, antialiased=False, alpha = 0.5)
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz,  color=rgba, zsort='average')
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.tick_params(axis='z', labelsize=10)
        ax.set_xlabel('N Primary', fontsize=10)
        ax.set_ylabel('N secondary', fontsize=10)
        ax.set_zlabel('Probability', fontsize=10)

    if flag[5]:
                    # Test Chan with GMMa and MC simulation
        t_et_TCA = 6e8
        x_p_MC_TCA = np.random.multivariate_normal(x_p_TCA, P_p_TCA, n_MC).T
        x_s_MC_TCA = np.random.multivariate_normal(x_s_TCA, P_s_TCA, n_MC).T
        GMM_p_TCA = []
        GMM_s_TCA = []
        PoC_out =np.zeros((n_splits_p, n_splits_s))
        #   GMM for the primary object
        for i in range(1,n_splits_p+1):
            GMM_p_TCA.append(GaussianMixture(n_components=i, tol = 1e-17).fit(x_p_MC_TCA.T))
        for i in range(1,n_splits_s+1):
            # GMM for the primary object
            GMM_s_TCA.append(GaussianMixture(n_components=i, tol = 1e-17).fit(x_s_MC_TCA.T))
        for i in range(0,n_splits_p):
            for j in range(0,n_splits_s):
                PoC_out[i, j] = PoC_ST(GMM_p_TCA[i],  GMM_s_TCA[j], R_sum,  PoC_order,  t_et_TCA, mu, flag_refinement)

        plt.figure()
        ax = plt.subplot(111, projection='3d')
        dx = np.ones(n_splits_p*n_splits_s)
        dy = np.ones(n_splits_p*n_splits_s)
        xpos = np.array([np.arange(0, n_splits_p, 1)]*n_splits_s).flatten()
        ypos = np.array([[i]*n_splits_p for i in range(n_splits_s)]).flatten()
        zpos = np.zeros(n_splits_p*n_splits_s)
        dz = PoC_out.flatten()
        X = np.arange(0, n_splits_p+1, 1)
        Y = np.arange(0, n_splits_p+1, 1)
        X, Y = np.meshgrid(X, Y)
        Z = np.ones((n_splits_p+1,n_splits_p+1))*PoC_true
        cmap = plt.cm.get_cmap('jet') # Get desired colormap - you can change this!
        max_height = np.max(dz)   # get range of colorbars so we can normalize
        min_height = np.min(dz)
            # scale each z to [0,1], and get their rgb values
        rgba = [cmap((k-min_height)/max_height) for k in dz] 
        surf = ax.plot_surface(X, Y, Z, color='green',
                    linewidth=0, antialiased=False, alpha = 0.5)
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz,  color=rgba, zsort='average')
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.tick_params(axis='z', labelsize=10)
        ax.set_xlabel('N Primary', fontsize=10)
        ax.set_ylabel('N secondary', fontsize=10)
        ax.set_zlabel('Probability', fontsize=10) 

    if flag[6]:
        t_et_TCA = 6e8
        GMM_pv = []
        w_p,v_p=np.linalg.eig(P_p_TCA)
        index_p = np.where(w_p == w_p.max())
        n=0

        for i in range(3,n_splits_p+1, 2):
            GMM_pv.append(GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[]))
            xi,Pi,wi = split_add(x_p_TCA, P_p_TCA, np.squeeze(v_p[:,index_p]), i)
            for j in range(i):
                GMM_pv[n].means_.append(xi[j, :])
                GMM_pv[n].covariances_.append(Pi[j, :, :])
                GMM_pv[n].weights_.append(wi[j])
            n +=1

        w_s,v_s=np.linalg.eig(P_s_TCA)
        index_s = np.where(w_s == w_s.max())
        GMM_sv =[]
        n = 0
        for i in range(3,n_splits_s+1, 2):
            GMM_sv.append(GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[]))
            xi,Pi,wi = split_add(x_s_TCA, P_s_TCA, np.squeeze(v_s[:,index_s]), i)
            for j in range(i):
                GMM_sv[n].means_.append(xi[j, :])
                GMM_sv[n].covariances_.append(Pi[j, :, :])
                GMM_sv[n].weights_.append(wi[j])
            n +=1

        n = 0
        n_el_p = len(GMM_pv)
        n_el_s = len(GMM_sv)
        PoC_outvv = np.zeros((n_el_p*n_el_s, ))
        for j in range(0,n_el_s):
            for i in range(0,n_el_p):
                PoC_outvv[n] = PoC_ST(GMM_pv[i],  GMM_sv[j], R_sum,  PoC_order,  t_et_TCA, mu, flag_refinement)
                n += 1

        plt.figure()
        ax = plt.subplot(111, projection='3d')
        dx = np.ones(n_el_p*n_el_s)
        dy = np.ones(n_el_p*n_el_p)
        xpos = np.array([np.arange(2, n_splits_p+1, 2)]*n_el_s).flatten()
        ypos = np.array([[i]*n_el_p for i in range(2, n_splits_s+1, 2)]).flatten()
        zpos = np.zeros(n_el_p*n_el_s)
        dz = PoC_outvv
        X = np.arange(2, 2+2*n_el_p)
        Y = np.arange(2, 2+2*n_el_s)
        X, Y = np.meshgrid(X, Y)
        Z = np.ones((2*n_el_s,2*n_el_p))*PoC_true
        cmap = plt.cm.get_cmap('jet') # Get desired colormap - you can change this!
        max_height = np.max(dz)   # get range of colorbars so we can normalize
        min_height = np.min(dz)
            # scale each z to [0,1], and get their rgb values
        rgba = [cmap((k-min_height)/max_height) for k in dz] 
        surf = ax.plot_surface(X, Y, Z, color='green',
                        linewidth=0, antialiased=False, alpha = 0.5)
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz,  color=rgba, zsort='average')
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.tick_params(axis='z', labelsize=10)
        ax.set_xlabel('N Primary', fontsize=10)
        ax.set_ylabel('N secondary', fontsize=10)
        ax.set_zlabel('Probability', fontsize=10)




plt.show()        



# elif flag == 2:
#     spice.furnsh('latest_leapseconds.tls')
#     n_split = 2
#     n_MC = 1000000
#     TCA_guess = 6e8
#     PoC_order = 30
#     mu = 398600.4415
#     flag_refinement = 1
#     R_sum = 15e-3

#     x_0_p = np.array([153951.4752631, 41874153.994752, 0,
#         3066.8746235984, -11.4110245828, 0])/1000 # [km/s]
#     x_0_s = np.array([153951.97327397, 41874156.744625, 2.752074327,
#     3066.8646233948, -0.0449985917, -11.3560272115])/1000 

#     P_0_p = np.array([[6494.0687699936000000000000, -376.2174293862200000000000, 0.0000000000000000000000,  0.0159835090297390000000 , -0.4942613967814000000000,   0.0000000000000000000000], 
#                         [-376.2174293862200000000000, 22.5686372188410000000000,0.0000000000000000000000 , -0.0009881643469155700000  , 0.0285754945542590000000 ,  0.0000000000000000000000],
#                         [ 0.0000000000000000000000,0.0000000000000000000000,1.2050405732107000000000 ,  0.0000000000000000000000 ,  0.0000000000000000000000,  -0.0000607087721201980000], 
#                         [-0.0009881643469155700000,0.0000000000000000000000  , 0.0000000443435166923780 , -0.0000012117808283202000,   0.0000000000000000000000, 0.0000607087721201980000],    
#                         [-0.4942613967814000000000,0.0285754945542590000000,0.0000000000000000000000 , -0.0000012117808283202000 ,  0.0000376229959779560000 ,  0.0000000000000000000000],    
#                         [ 0.0000000000000000000000,0.0000000000000000000000,-0.0000607087721201980000 ,  0.0000000000000000000000  , 0.0000000000000000000000  , 0.0000000033903879281480]])/1000**2


#     P_0_s = np.array([[  6539.7159476886000000000000,-354.4827622708000000000000 ,  -24.2158339268500000000000 ,  0.0161364407930240000000 , -0.4975675313434400000000 ,-0.0000667805505728000000], 
#                         [  -354.4827622708000000000000,19.9876841613490000000000,1.3128407317001000000000 , -0.0009366299987872100000  , 0.0269117100210450000000 ,  0.0000038371255372188000],
#                         [-24.2158339268500000000000,1.3128407317001000000000,1.2674955658292000000000  ,-0.0000599870075272560000 ,  0.0018427685478909000000 , -0.0000602391913489130000], 
#                         [0.0161364407930240000000,-0.0009366299987872100000,-0.0000599870075272560000 ,  0.0000000447836341841990  ,-0.0000012229871227354000 , -0.0000000001697954617528],    
#                         [-0.4975675313434400000000,0.0269117100210450000000,0.0018427685478909000000 , -0.0000012229871227354000   ,0.0000378619099092410000  , 0.0000000050463685349191],    
#                         [ -0.0000667805505728000000,0.0000038371255372188000,-0.0000602391913489130000 , -0.0000000001697954617528  , 0.0000000050463685349191 ,  0.0000000034465981691675]])/1000**2

#     # Test Chan at TCA
#     x_0_p_MC = np.random.multivariate_normal(x_0_p, P_0_s, n_MC).T
#     GMM_p = GaussianMixture(n_components=n_split, tol = 1e-17).fit(x_0_p_MC.T)

#     # GMM for the primary object
#     x_0_s_MC = np.random.multivariate_normal(x_0_s, P_0_s, n_MC).T
#     GMM_s = GaussianMixture(n_components=n_split, tol = 1e-17).fit(x_0_s_MC.T)
#     GMM_TCA = namedtuple("GMM_TCA" , "means_ covariances_ weights_ TCA_et_ TCA_UTC_")
#     # GMM_p = GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[])
#     # GMM_s = GMM_TCA(means_= [], covariances_ = [], weights_ = [], TCA_et_=[], TCA_UTC_=[])
#     # GMM_p.means_.append(x_0_p)
#     # GMM_s.means_.append(x_0_s)
#     # GMM_p.covariances_.append(P_0_p)
#     # GMM_s.covariances_.append(P_0_s)
#     # GMM_p.weights_.append(1)
#     # GMM_s.weights_.append(1)
#     # Output PoC
#     PoC_out= PoC_ST(GMM_p, GMM_s, R_sum,  PoC_order, TCA_guess, mu, flag_refinement)
    
#     print('c')
