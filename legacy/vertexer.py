""" 
vertexer copied from 
 https://codereview.stackexchange.com/questions/278487/solving-the-tdoa-multilateration-problem-in-3-dimensions 
users 10Gev and Reinderien

Based on multilateration navigation approach: https://en.wikipedia.org/wiki/Trilateration

"""
from dataclasses import dataclass
import scipy
from numpy.random import default_rng
from scipy.optimize import minimize
from scipy import constants
import numpy as np


@dataclass
class Vertexer:
    roc: np.ndarray
    c: float

    def lin_est(self, times: np.ndarray) -> np.ndarray:
        system = np.hstack((self.roc, np.atleast_2d(times).T))
        LHS = np.tile(system, (len(system), 1))
        RHS = np.repeat(system, len(system), 0)

        A = 2*(LHS - RHS)
        b = np.square(LHS) - np.square(RHS)
        c2 = self.c**2
        A[:, -1] *= c2
        b[:, -1] *= c2
        b = b.sum(axis=1)

        return np.linalg.lstsq(A, b, rcond=None)[0]

    def objective(self, var: np.ndarray, times: np.ndarray) -> float:
        norms = np.linalg.norm(var[:3] - self.roc, axis=1)
        addends = norms + self.c*(var[3] - times)
        chi2 = addends.sum()
        return chi2

    def find(self, init: np.ndarray, times: np.ndarray) -> np.ndarray:
        res = minimize(fun=self.objective, x0=init, args=times,
                       method='COBYLA', options={'maxiter': 1e5})
        return res.x.flatten()



def main() -> None:

    c = constants.c
    tx_XYZ = np.array([[  4057484.7830325 ,  14568429.37210906,  21603872.89041781],
       [-23901103.55812721,  10626663.89092992,  -3352745.77814737],
       [-21554012.90592577, -10265229.75241027,  11201751.05871473],
       [-20160012.44460722,   2210182.88030742,  16855818.90235173],
       [ -5621253.99235301,  18487399.50257921,  17945673.68993121],
       [ -2056471.42532162,  26290410.10973211,   -122075.95558358],
       [-14444946.68297005,  12045728.68172588,  18484414.43087659]])
    rx_XYZ = np.array([-3766938.45081918,  3211854.8708333 ,  4008305.67272309])
    obs_delays_sec = np.array([0.07457209, 0.0756646 , 0.07821047, 0.0695539 , 0.0692522 , 0.07841279, 0.06684736])  # typical ionosphere error
    # obs_delays_sec = np.linalg.norm(rx_XYZ - tx_XYZ, axis=1) / c  # Zero ionosphere error
    my_vertexer = Vertexer(-tx_XYZ, c)

    print('Transmitter coordinates -')
    print('Actual:', rx_XYZ)

    init = my_vertexer.lin_est(obs_delays_sec)
    print('Initial-estimate:', init)

    x = my_vertexer.find(init, obs_delays_sec)
    print('Estimated: ', x)

    print('Linear err: %1.1e m\nMinimized err: %1.1e m\n' %\
         (np.linalg.norm(init[:3] - rx_XYZ), np.linalg.norm(x[:3] - rx_XYZ)))



if __name__ == '__main__':
    main()
