import numpy as np
import math
import matplotlib.pyplot as plt


class Material:
    def __init__(self, concentration, temperator, Ea):
        kb = 1.381e-23 #[J/K]
        Na = 6.022e23 #[1/mol]
        R = 8.314 #[J/mol*K]
        self.concentration = concentration
        self.k = (kb/Na)*temperator*np.exp(-(Ea/(R*temperator)))

def timestep(A, B):
    dt = 0.01/(max(A.concentration, B.concentration)*max(A.k, B.k))
    A_con = dt*(-A.k*A.concentration+B.k*B.concentration) + A.concentration
    B_con = dt*(-B.k*B.concentration+A.k*A.concentration) + B.concentration
    A.concentration = A_con
    B.concentration = B_con
    return dt


if __name__ == '__main__':
    kcal_to_J = 1000*4.184

    Ea = 25*kcal_to_J #[J/mol]
    DeltaG = 2*kcal_to_J #[J/mol]
    tmp = [200, 300, 400]
    A0 = [1, 0.5, 0, 1]
    B0 = [0, 0.5, 1, 0.5]

    def plotABtime(A0, B0, tmp, Ea, DeltaG):
        A = Material(A0, tmp, Ea)
        B = Material(B0, tmp, Ea+DeltaG)
        time = [0]
        A_lst = [A.concentration]
        B_lst = [B.concentration]
        for i in range(10000):
            time.append(timestep(A, B) + time[i])
            A_lst.append(A.concentration)
            B_lst.append(B.concentration)
        plt.figure
        plt.plot(time, A_lst, label=f'A - {tmp}K, A0 - {A0}, B0 - {B0}')
        plt.plot(time, B_lst, label=f'B - {tmp}K, A0 - {A0}, B0 - {B0}')

    for t in tmp:
        for i in range(len(A0)):
            plotABtime(A0[i], B0[i], t, Ea, DeltaG)



    # Customize axis labels
    plt.xlabel('Concentration')
    plt.ylabel('Time')

    plt.xscale('log')

    # Add a legend
    plt.legend()

    # Add a title
    plt.title('A to B to eq')

    # Display the plot
    plt.show()
