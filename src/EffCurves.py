import argparse 
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

"""_summary_

We used this piece of code to calculate efficency curve for each PMT and find a good working point.

The parser is a folder name, with files called PMT(number of PMT).txt. 

In this file we will assume a familiar structure, which is the one the owner of this file used in the first place
"""


def eff_curve(data, A = 1 , V_ind = 0, T_ind = 1, C_ind = 2, PMTu_ind = 3, PMTd_ind = 4):
    
    """_summary_
    
    The function estimates efficency as ratio between triples and couples, corrects with acceptance and uses single count to remove fake couples.
    
    Error is estimated from binomial stats, the users will recognize the variance of a binomial distribution
    
    Args:
        data (ndarray) : is a data numpy array,here assumed by the following structure : HV/Triples/Couples/Singles of upper PMT/Singles of bottom PMT  
        A (float, optional) : Acceptance, is a scale geometric factor for efficency. Defaults to 1.
        V_ind (int, optional) : HV index. Defaults to 0.
        T_ind (int, optional) : Triple index. Defaults to 1.
        C_ind (int, optional) : Couples index. Defaults to 2.
        PMTu_ind (int, optional) : Upper PMT singles index. Defaults to 3.
        PMTd_ind (int, optional) : Bottom PMT singles index. Defaults to 4.

    Returns:
        eff(ndarray) : efficency estimated as explained above.
        eff_err(ndarray) : efficency error estimated as explained above.
        C_fake(ndarray) : fake couples estimated as product of singles divided by 100 s multiplied by 20 ns
    """
    
    HV = data[:, V_ind] * 1e-3
    
    T = data[:, T_ind]
    
    C = data[:, C_ind]
    
    S_07 = data[:, PMTu_ind]
    
    S_02 = data[:, PMTd_ind]
    
    C_fake = S_07 * S_02 * 8e-10 
    
    eff = T/(A*(C-C_fake))
    
    eff_err = np.sqrt(eff * (1-eff)/ (C-C_fake))
    
    return HV, eff, eff_err, C_fake

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "This program accepts a path to a folder, it estimates efficency and associated error, and plots efficency as a function of HV")

    parser.add_argument("folder", type = str, help = "Path.txt")

    args = parser.parse_args()

    folder = args.folder

    files = glob.glob(os.path.join(folder, "PMT*.txt"))

    A = {"PMT01": 1, "PMT02" : 0.4, "PMT04" : 1 , "PMT07" : 0.5 , "PMT08" : 1, "PMT09" : 1, "PMT10" : 1, "PMT11" : 1,"PMTOR" : 1}
    
    fig, ax = plt.subplots()
    
    ax.set_xlabel('HV [kV]')
    
    ax.set_ylabel('Efficiency [Pure]')
    
    ax.grid(True)
    
    labels = []

    for f in files:
        title = os.path.splitext(os.path.basename(f))[0]
        
        HV, eff, err, C_fake = eff_curve(np.loadtxt(f), A[title])
        
        ax.errorbar(HV, eff, err, fmt = '.' , capsize = 3)
        
        labels += [title]
        
    ax.legend(labels)
    
    folder_name = os.path.basename(os.path.normpath(folder))
    
    ax.set_title(f'Efficiency Curves - {folder_name}')
    
    plt.savefig(f'Efficiency Curves - {folder_name}')
    
    plt.show()
    