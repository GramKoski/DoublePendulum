import numpy as np
import sympy as smp
from scipy.integrate import solve_ivp
from timeit import default_timer as timer


assemble = lambda S, w: np.hstack((S,w.flatten()))
disassemble = lambda state: (state[:4], state[4:].reshape(4, 4))



def getLyap(dSdt, jacobian_f, t, inState, args):
    start = timer()
    g, m1, m2, L1, L2 = args
    S = inState

    def expandedFunc(t, state, *args):
        S, w = disassemble(state)
        t1, t1d, t2, t2d = S
        S_d = dSdt(t, S, *args)
        w_d = jacobian_f(t, *args, t1, t1d, t2, t2d) @ w
        return assemble(S_d, w_d)
    
    dt = 1
    iters = 40
    lyaps = []
    w = np.identity(4)

    for _ in range(iters):
        sol = solve_ivp(expandedFunc, [0, dt], assemble(S,w), t_eval=[dt], args=(g, m1, m2, L1, L2), max_step=dt)
        S,w = disassemble(sol.y.flatten())
        w,r = np.linalg.qr(w)
        lyaps.append(np.log(abs(r.diagonal()))/dt)

    transient_steps = 0
    maxlyap = max(np.average(lyaps[transient_steps:],axis=0))
    end = timer()
    print('Lambda = ' + str(maxlyap) + ', time_elapsed = ' + str(end - start))
    return maxlyap

'''
def getLyap(dSdt, inState, *args):
    start = timer()

    dt = 0.01
    totalTime = 40
    checkpoint = 0.2
    g, m1, m2, L1, L2 = args
    S = inState

    theta1, theta2, the1_d, the2_d = inState 
    perturbed = [theta1+0.001, the1_d, theta2, the2_d]
    initialD = np.linalg.norm(np.subtract(perturbed, S))
    exponents = []
    

    for i in range(int(totalTime/checkpoint)):
        S = solve_ivp(dSdt, t_span = [0, checkpoint], y0 = S, t_eval = [checkpoint], args=(g, m1, m2, L1, L2), max_step=checkpoint).y.flatten()
        perturbed = solve_ivp(dSdt, [0, checkpoint], perturbed, t_eval = [checkpoint], args=(g, m1, m2, L1, L2), max_step=checkpoint).y.flatten()
        finalD = np.linalg.norm(np.subtract(perturbed, S))
        exponents.append(np.log(abs(finalD/initialD)))
        perturbed = S + (initialD/finalD)*(np.subtract(perturbed,S))
        initialD = np.linalg.norm(np.subtract(perturbed,S))

    end = timer()
    lyap = np.sum(exponents)/totalTime
    print('Lambda = ' + str(lyap) + ', time_elapsed = ' + str(end - start))
    return lyap
'''











