from glob import glob
import mat73
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import zscore
from scipy.sparse import csr_matrix
from scipy.ndimage import gaussian_filter1d
from natsort import natsorted
import prediction, datasets
import imp, sys
from natsort import natsorted

def summary_stats(root, mouse, day, syncope_times):
    syncope, fnames = syncope_times
    dat = np.load(root + f'proc/dat_{mouse}_D{day}.npy', allow_pickle=True).item()
    proc = np.load(root + f'proc/procnn_{mouse}_D{day}.npy', allow_pickle=True).item()
    inds = np.nonzero(np.array([mouse in f and ('D'+str(day) in f or 'Day'+str(day) in f) for f in fnames]))[0]
    syncope_on = syncope['thresh'][:,inds[-1]]
    syncope_off = syncope['thresh2'][:,inds[-1]]

    S = dat['neurons']['S']
    S_mean = S.mean(axis=1).copy()[:,np.newaxis]
    S_std = S.std(axis=1).copy()[:,np.newaxis]
    print('S = ', S.shape)
    sareas, area_names = dat['neurons']['sareas'], dat['neurons']['area_names']

    S_pred_test = proc['S_pred_test']
    S_pred_test_neigh = proc['S_pred_test_neigh']
    S_test = proc['S_test']    
    firing_rate = (proc['S_val'] * 10).mean(axis=1)
    firing_variance = proc['S_val'].var(axis=1)

    n_neurons = S_test.shape[0]
    S_test = S_test.reshape(n_neurons, -1, 610)
    S_pred_test = S_pred_test.reshape(n_neurons, -1, 610)
    S_pred_test_neigh = S_pred_test_neigh.reshape(n_neurons, -1, 610)
    n_stims = S_test.shape[1]

    #print(len(proc['laser_durations']), S_test.shape[1], proc['laser_hz'], proc['laser_durations'])
    # take 20Hz 30s trial
    #ind = np.nonzero(np.logical_and(proc['laser_durations'] < 50, proc['laser_hz'] == 20))[0][0]
    if proc['laser_hz'][-1] == 20:
        inds = np.arange(0, n_stims)[-4:]
    else:
        inds = np.arange(0, n_stims)[1:-1]
    S_test = S_test[:, inds]
    S_pred_test = S_pred_test[:, inds]
    S_pred_test_neigh = S_pred_test_neigh[:, inds]

    from scipy.stats import ttest_rel
    resp_all = np.zeros((0,))

    l0 = 10
    residual = S_test.copy() - S_pred_test_neigh.copy()
    dtlas = 8
    laser_on = residual[:, :, l0 : l0+dtlas].mean(axis=-1)
    laser_pre = residual[:, :, l0-1-dtlas : l0-1].mean(axis=-1)
    #resp = laser_on - laser_pre

    p = np.array([ttest_rel(laser_on[i], laser_pre[i], alternative='greater').pvalue 
          for i in range(len(laser_on))])
    p2 = np.array([ttest_rel(laser_on[i], laser_pre[i], alternative='less').pvalue 
          for i in range(len(laser_on))])

    responsive = np.logical_and(laser_on.mean(axis=-1) > 0, p < 0.05)

    laser_fr = S_test[:, :, l0 : l0+dtlas].mean(axis=-1)

    #responsive = np.logical_and(responsive, 
    #                            laser_fr.mean(axis=-1)*10 > firing_rate + 1)

    print(f' \t {responsive.mean(): 0.3f} \t {(p < 0.05).mean(): 0.3f} \t {(p2 < 0.05).mean(): 0.3f}') #{(resp > 0.5).mean(): 0.3f}, {(resp < -0.5).mean(): 0.3f}
    #resp_all = np.append(resp_all, p2[isort])


    files = natsorted(glob(root + f'Data/*{mouse}*_D{day}*_imec*'))
    files = files if len(files) > 0 else natsorted(glob(root + f'Data/*{mouse}*_Day{day}_g0_imec*'))
    print(files)
    dats = [mat73.loadmat(f) for f in files]

    st = np.zeros(0, 'float64')
    clu = np.zeros(0, 'uint32')
    nc0 = 0
    for i in range(len(dats)):
        cluid = dat['neurons']['clu_id'][dat['neurons']['probe_id']==i]
        good_spks = np.isin(dats[i]['sp']['clu'], cluid)
        clui = dats[i]['sp']['clu'][good_spks]
        #clui = np.unique(clui, return_inverse=True)[1]
        st = np.append(st, dats[i]['sp']['st'][good_spks], axis=0)
        clu = np.append(clu, dats[i]['sp']['clu'][good_spks] + nc0, axis=0)
        if i==0:
            nc0 = dats[i]['sp']['clu'].max() + 1

    clu_id = dat['neurons']['clu_id']
    probe_id = dat['neurons']['probe_id']
    clu_id[probe_id==1] += nc0    

    digArray = dats[0]['digArray']
    laser_times = np.unique(np.nonzero((np.diff(digArray[2,:])>0))[0])
    laser_times = np.hstack((laser_times[0], laser_times[1:][np.diff(laser_times) > 5 * 12500/10])).astype('float64')
    laser_times /= 12500
    laser_control = laser_times[0]
    laser_times = laser_times[inds]
    print(laser_times)

    dt = 0.1
    t0 = -30.0
    t1 = 90.0
    tpts = np.arange(t0, t1+dt, dt)
    npts = len(tpts)-1
    raster_summary = np.nan * np.ones((n_neurons, npts), 'int')
    lt = laser_times[-2] # 20Hz 30s stim laser on
    xrange = lt + np.array([t0, t1])
    inrange = np.logical_and(st > xrange[0], st < xrange[1])
    n_clu = clu.max()+1
    n_neurons = S_test.shape[0]
    tfs_all = -1 * np.ones((n_neurons,))
    sts = st[inrange].copy() - lt
    clus = clu[inrange]

    for c in range(n_clu):
        cs = np.nonzero(clu_id == c)[0]
        if len(cs) > 0:
            iclu = clus == c
            stb = np.histogram(sts[iclu], tpts)[0]
            raster_summary[cs[0]] = stb
    
    dt = 0.001
    t0 = -2.0
    synct = syncope_off[0] - syncope_on[0]
    t1 = synct 
    nl = int(t0/dt)*-1
    tpts = np.arange(t0, t1+dt, dt)
    npts = len(tpts)-1

    loffset = [syncope_on[0], # timeoff, inactive
               syncope_on[0], # timeoff, inactive control
               0.0, #timeoff, inactive from start of laser
               syncope_off[0], # latency to response after syncope
               max(30.0, syncope_off.max())] # latency to response after laser
    tend = [t1, t1, 30.0, 6.0, 6.0]
    print(loffset)

    for ll, laser_time in enumerate([laser_times[-2], laser_control, laser_times[-2], laser_times[-2], laser_times[-2]]):
        lt = laser_time + loffset[ll]
        xrange = lt + np.array([t0, t1])
        inrange = np.logical_and(st > xrange[0], st < xrange[1])
        n_clu = clu.max()+1
        n_neurons = S_test.shape[0]
        tfs_all = -1 * np.ones((n_neurons,))
        sts = st[inrange].copy() - lt
        clus = clu[inrange]

        plt.figure(figsize=(8,8))
        plt.scatter(sts, clus, s=1)
        plt.plot([0,0], [0,clus.max()])
        plt.plot([t1-0.5,t1-0.5], [0,clus.max()])

        raster = np.nan * np.ones((n_neurons, 2500), 'int')
        latency = np.nan * np.ones((n_neurons,), 'float32')
        if ll==0:
            timeoff = np.zeros((n_neurons,), 'float32')
        for c in range(n_clu):
            cs = np.nonzero(clu_id == c)[0]
            if len(cs) > 0:
                iclu = clus == c
                stb = np.histogram(sts[iclu], tpts)[0][:-1]

                fr = firing_rate[cs[0]]
                N = 4
                #ppoisson = np.ones_like(stb)
                nsp0 = 0
                dtf = min(int(4/dt), int(-np.log(0.01)/fr / dt))
                raster[cs[0]] = stb[:2500]
                
                if ll < 3:
                    nolat = True
                    for t in range(nl, len(stb)-dtf, dtf):
                        nsp = stb[t : t+dtf].sum()
                        if ll==0:
                            timeoff[cs[0]] += (nsp==0) * dtf * dt
                        #print(dtf)
                        if nsp < 1 and nolat:
                            latency[cs[0]] = (t - nl) * dt
                            plt.scatter(latency[cs[0]], c, color='k', s=5)
                            nolat = False
                else:
                    for t in range(nl, len(stb)):
                        nsp = stb[nl : t+1].sum()
                        if nsp > 0 and nsp0!=nsp:
                            n = np.arange(0, nsp-1)
                            try:
                                ppoisson = 1 - np.array([(N * fr * tpts[t+1-nl])**m * 
                                                        np.exp(- N * fr * tpts[t+1-nl]) / np.math.factorial(m) 
                                                         for m in n]).sum()
                            except:
                                ppoisson = 0
                            nsp0 = nsp
                            if ppoisson < 1e-6:
                                latency[cs[0]] = (t - nl + 1) * dt
                                plt.scatter(latency[cs[0]], c, color='k', s=2)
                                break
                    
        if ll==0:
            timeoff /= synct

            syncon = int(-t0*10 + syncope_on[0]*10)
            syncoff = int(-t0*10 + syncope_off[0]*10)
            res_sync = residual[:, 3, syncon : syncoff]
            #res_rand = S_val - S_val_pred_neigh[:, 3, syncon : syncoff]

            rs = res_sync.mean(axis=1)
            rs[rs >= 0] = np.nan
            rs0 = residual[:,:,int(-t0*10)-10:int(-t0*10)-1].mean(axis=(-2,-1))

            suppressed = np.logical_and(rs < rs0, latency < 0.25)
            print(suppressed.mean(), (latency<0.25).mean())

            latencyS = latency.copy()
            rasterS = raster.copy()
        elif ll==1:
            latencyC = latency.copy()
            rasterC = raster.copy()
        elif ll==2:
            latencyS0 = latency.copy() - syncope_on[0]
            rasterS0 = raster.copy()
        elif ll==3:
            latencyR = latency.copy()
            rasterR = raster.copy()
        else:
            latencyO = latency.copy()
            rasterO = raster.copy()
        
            

        plt.show()

    dt = 0.0025
    t0 = -0.1
    t1 = 0.5
    nl = int(t0/dt)*-1
    tpts = np.arange(t0, t1+dt, dt)
    npts = len(tpts)-1        

    sts = np.zeros((0,), 'float64')
    clus = np.zeros((0,), 'uint32')
    n_clu = clu.max()+1
    frs = np.zeros((n_clu, 4), 'float32')
    for i in range(0,4):

        lt = laser_times[i]
        xrange = lt + np.array([t0, t1])
        inrange = np.logical_and(st > xrange[0], st < xrange[1])
        tfs_all = -1 * np.ones((n_neurons,))

        sts = np.append(sts, st[inrange] - lt, axis=0)
        clus = np.append(clus, clu[inrange], axis=0)
        
        lt = laser_times[i]
        xrange = lt + np.array([0, 2.0])
        inrange = np.logical_and(st > xrange[0], st < xrange[1])
        clu0 = clu[inrange]
        nspks = np.histogram(clu0, np.arange(0, n_clu+1)-0.5)[0]
        frs[:,i] = nspks / 2.0

    plt.figure(figsize=(8,8))
    plt.scatter(sts, clus, s=1)

    raster = np.nan * np.ones((n_neurons, npts-1), 'float32')
    latency = np.nan * np.ones((n_neurons,), 'float32')
    laser_fr = np.nan * np.ones((n_neurons,4), 'float32')
    for c in range(n_clu):
        cs = np.nonzero(clu_id == c)[0]
        if len(cs) > 0:
            laser_fr[cs[0]] = frs[c]
        if len(cs) > 0 and responsive[cs[0]]:
            iclu = clus == c
            stb = np.histogram(sts[iclu], tpts)[0][:-1]
            raster[cs[0]] = stb
            laser_fr[cs[0]] = frs[c]
            fr = firing_rate[cs[0]]
            N = 4
            #ppoisson = np.ones_like(stb)
            nsp0 = 0
            for t in range(nl, len(stb)):
                nsp = stb[nl : t+1].sum()
                if nsp > 0 and nsp0!=nsp:
                    n = np.arange(0, nsp-1)
                    ppoisson = 1 - np.array([(N * fr * tpts[t+1-nl])**m * 
                                            np.exp(- N * fr * tpts[t+1-nl]) / np.math.factorial(m) for m in n]).sum()
                    nsp0 = nsp
                    if ppoisson < 1e-6:
                        latency[cs[0]] = (t - nl + 1) * dt
                        plt.scatter(latency[cs[0]], c, color='k', s=2)
                        break

    plt.title(f'{mouse} D{day}')
    plt.show()

    print(latency[latency>0].mean())
    latencyL = latency.copy()
    rasterL = raster.copy()
    
    brain_areas = proc['area_names'][proc['sareas']]
    
    return (firing_rate, firing_variance, raster_summary, suppressed, laser_fr, 
            latencyS, rasterS, latencyO, rasterO, latencyC, rasterC, 
            latencyS0, rasterS0, latencyR, rasterR, timeoff, 
            latencyL, rasterL, brain_areas)
