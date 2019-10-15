import numpy as np
import lmfit as lm
import profiles

def calcmedquartnine(array):
    index = np.argsort(array,axis=0)
    median = np.median(array)
    sixlow = np.percentile(array, 16)
    sixhigh = np.percentile(array, 84)
    ninelow = np.percentile(array, 2.5)
    ninehigh = np.percentile(array, 97.5)
    nineninelow = np.percentile(array, 0.15)
    nineninehigh = np.percentile(array, 99.85)

    return median, sixlow, sixhigh, ninelow, ninehigh,\
        nineninehigh, nineninelow



def residual(params, x, data, eps_data):
    M1 = params['M1'].value
    M2 = params['M2'].value
    M3 = params['M3'].value
    a1 = params['a1'].value
    a2 = params['a2'].value
    a3 = params['a3'].value

    model = profiles.threeplumsurf(x,M1,M2,M3,a1,a2,a3)

    return (data-model)/eps_data


def tracerfit(p0in,p0in_min,p0in_max,\
              rbin_phot,surfden,surfdenerr):
    params = lm.Parameters()
    params.add('M1', value=p0in[0],min=p0in_min[0],max=p0in_max[0])
    params.add('M2', value=p0in[1],min=p0in_min[1],max=p0in_max[1])
    params.add('M3', value=p0in[2],min=p0in_min[2],max=p0in_max[2])
    params.add('a1', value=p0in[3],min=p0in_min[3],max=p0in_max[3])
    params.add('a2', value=p0in[4],min=p0in_min[4],max=p0in_max[4])
    params.add('a3', value=p0in[5],min=p0in_min[5],max=p0in_max[5])

    out = lm.minimize(residual, params, \
                      args=(rbin_phot, surfden, surfdenerr),\
                      maxfev=500)

    pfits = np.zeros(6)
    pfits[0] = out.params['M1'].value
    pfits[1] = out.params['M2'].value
    pfits[2] = out.params['M3'].value
    pfits[3] = out.params['a1'].value
    pfits[4] = out.params['a2'].value
    pfits[5] = out.params['a3'].value

    return pfits

def v4func(x,A,B,C):
    return (A)*x**B + (C)

def v4residual(params, x, data, eps_data):
    A = params['A'].value
    B = params['B'].value
    C = params['C'].value

    model = v4func(x,A,B,C)

    return (data-model)/eps_data

def v4fit(p0in,p0in_min,p0in_max,\
          rbin,v4,v4err):
    params = lm.Parameters()
    params.add('A', value=p0in[0],min=p0in_min[0],max=p0in_max[0])
    params.add('B', value=p0in[1],min=p0in_min[1],max=p0in_max[1])
    params.add('C', value=p0in[2],min=p0in_min[2],max=p0in_max[2])

    out = lm.minimize(v4residual, params, \
                      args=(rbin, v4, v4err))

    pfits = np.zeros(3)
    pfits[0] = out.params['A'].value
    pfits[1] = out.params['B'].value
    pfits[2] = out.params['C'].value

    return pfits


def residual_powline(params, x, data, eps_data):
    A = params['A'].value
    B = params['B'].value
    C = params['C'].value

    model = A*x**B+C

    return (data-model)/eps_data


def fit_powline(rbin_tmp,vlos4med,vlos4err,Rhalf):
    max_r = np.max(rbin_tmp)
    j=0L
    if (max_r > Rhalf):
        while(rbin_tmp[j] < Rhalf):
            j=j+1
    else:
        j = len(rbin_tmp)-4 # just go back 4 bins if data ends at rhalf
    ruse_t = rbin_tmp[j:]/max_r   #normalised by rdata
    vuse_t = vlos4med[j:]
    verr_t = vlos4err[j:]

    #Upsample:
    ruse = np.logspace(np.log10(np.min(ruse_t)),np.log10(np.max(ruse_t)),\
                       100)
    vuse = np.interp(ruse,ruse_t,vuse_t)
    verr = np.interp(ruse,ruse_t,verr_t)

    params_powline = lm.Parameters()
    params_powline.add('A', value=(np.max(vlos4med)),\
               min=(np.min(vlos4med)/100.0),\
               max=(np.max(vlos4med)*100.0))
    params_powline.add('B', value=0.0,\
               min=-2.0,\
               max=2.0)
    params_powline.add('C', value=0.0,vary=False)

    out = lm.minimize(residual_powline, params_powline, \
                      args=(ruse, vuse, verr))
    pfits_powline = np.zeros(3)
    pfits_powline[0] = out.params['A'].value
    pfits_powline[1] = out.params['B'].value
    pfits_powline[2] = out.params['C'].value

    #print "Power law fit", pfits_powline
    return pfits_powline


def tvl4func(rint,rbin_tmp,vlos4med,A,B,C,rout,gamout):
    tvl4 = np.interp(rint,rbin_tmp,vlos4med,left=0,right=0)

    for i in range(len(tvl4)):
        if (rint[i] > np.max(rbin_tmp)):
            tvl4[i] = A*(rint[i]/np.max(rbin_tmp))**B+C
        if (rint[i] > rout):
            tvl4[i] = (A*(rout/np.max(rbin_tmp))**B+C)*\
                (rint[i]/rout)**(-gamout)
    return tvl4
