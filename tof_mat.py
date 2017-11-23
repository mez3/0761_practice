from numpy import *
import scipy.io
from time import time
from multiprocessing import Pool
import scipy.signal
import matplotlib.pyplot as mpl

MED = 50
TOFCF_ID_IDX = 0
delay_cf = 1000
shft = 3000
etalon = scipy.io.loadmat('detector_filtered_field.mat')
d = scipy.io.loadmat('det_struct_2.mat')
det_struct = array(d['det_struct'])
n = det_struct.shape[1]
data_et = array(etalon['detector_filtered_field'])
s = array(det_struct['detect'])

def tof_mat(start_k_amp, stop_k_amp, p_filt):
    pool = Pool(8)
    args = list(zip(range(2,n),[start_k_amp]*(n-1),[stop_k_amp]*(n-1), [p_filt]*(n-1)))
    det_CF_m = list(map(main_1_calc,args))
    #for j in range(2,n):
        #det_CF_m[j-1] = main_1_calc(j,start_k_amp,stop_k_amp)
    print(visualize(det_CF_m))

def main_1_calc(args):
    j = args[0]
    start_k_amp = args[1]
    stop_k_amp = args[2]
    p_filt=args[3]
    len_mas = stop_k_amp - start_k_amp +1 
    data = s[0][j]
    data2 = signal_cut(data)
    id_cf_mas = zeros((int8(len_mas),1))
    u=1
    for q in range(start_k_amp, stop_k_amp-1):
        id_cf_mas[u-1] = tof_cf(data2,delay_cf,TOFCF_ID_IDX,q)
        u=u+1
    mpl.show()
    if p_filt > stop_k_amp - start_k_amp:
        if (stop_k_amp -start_k_amp) % 2 == 0:
            p_filt = stop_k_amp - start_k_amp - 3
        else:
            p_filt = stop_k_amp - start_k_amp - 2
    id_cf_mas = scipy.signal.medfilt(id_cf_mas, p_filt)
    id_cf = median(id_cf_mas)
    return id_cf

def signal_cut(s):
    shift = 1000
    s_shift = s[shift:-1]
    s_fill = s[shift] * ones((shift-1,1))
    return vstack((s_fill, s_shift))

def constant_fraction(y,delay):
    f = 0.5
    #EndVal=y.max()
    #print(EndVal)
    offset = y[delay-1]*ones((delay,1))
    offset2= y[-1]*ones((delay,1))
    f1 = vstack((offset, y))
    f2 = f * vstack((y, offset2))
    f_res = f1 - f2
    return f_res


def tof_cf(s,delay,ID_IDX,k_amp): 

    Ythreshold = k_amp * max(s) / 100
    
    y = s - Ythreshold
  
    f_res = constant_fraction(y,delay)

    # % Zerro-crossing %
    t1 = f_res.shape[0]
    k = 0
    idn = []
    for i in range(1 + delay, t1):
        if f_res[i] * f_res[i - 1] < 0:
            idn.append(i - delay)
            k = k + 1
   # %  Making id_cf 
    k = k - 1
    if any(idn):
        if ID_IDX > k:
            if idn[k-1] > 1:
                id_cf = idn[k-1]
            else:
                id_cf = delay
        else:
            id_cf = idn[ID_IDX]
    else:
        id_cf = y.shape[0]
            
    if id_cf >= f_res.shape[0]:
        id_cf=y.shape[0]
    else:
        id_cf=id_cf-1
    #id_cf = 0
    #for i in range(0, len(idn)):
        #if idn[i] + delay < sum(nonzero(f_res == f_res.max())[0]) / len(nonzero(f_res == f_res.max())[0]) and idn[i] + delay > sum(nonzero(f_res == f_res.min())[0])/len(nonzero(f_res == f_res.min())):
            #id_cf = idn[i]
    return id_cf

def visualize(det_CF_m):
    RANGE=450
    m_y2 = median(det_CF_m[2:18])
    er_y2_1 = nonzero(det_CF_m < (m_y2-RANGE))
    er_y2_2 = nonzero(det_CF_m > (m_y2+RANGE))
    errors2 = hstack((er_y2_1, er_y2_2))
    n_er = errors2.shape[1]
    k_er = n_er * 100 / len(det_CF_m)
    return k_er

if __name__ == "__main__":
    tostart = time()
    args = [35,40,11]
    tof_mat(*args)
    print('%f m %f s'% ((time() - tostart)//60, (time() - tostart)%60))
