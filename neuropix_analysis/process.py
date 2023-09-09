import numpy as np
from rastermap import Rastermap
import mat73
from rastermap.utils import bin1d
from scipy.ndimage import gaussian_filter1d, maximum_filter1d, minimum_filter1d, median_filter
from scipy.stats import skew
import torch
from glob import glob
import datasets, prediction, summary

mice = ['18N', '22R', '32N', '47L', '57L', '59L', '74L']
days = [[1], [1], [1,2], [1,2], [1,2], [1,2], [1,2]]

def process_dataset(root, mouse, day, device=torch.device("cuda")):
    
    dat = datasets.load_neurons_behaviors(root, mouse, day)

    S = dat['neurons']['S']
    S_mean = S.mean(axis=1).copy()[:,np.newaxis]
    S_std = S.std(axis=1).copy()[:,np.newaxis]
        
    print('S = ', S.shape)
    sareas, area_names = dat['neurons']['sareas'], dat['neurons']['area_names']

    # sort neurons to visualize data
    bin_size = 5 # bin in neurons
    rm = Rastermap(time_lag_window=10, 
                   n_clusters=50 if S.shape[0] > 300 else 25, 
                   locality=0.75, 
                   n_PCs=200 if S.shape[0] > 300 else 100,
                   normalize=False, 
                   mean_time=False,
                   bin_size=bin_size,
                   verbose=False).fit((S - S_mean) / S_std,  # use z-scored responses to get PCs
                                       )
    print('n_clusters = ', rm.embedding_clust.max() + 1)
    # PCs of data
    U = rm.Usv.copy()
    U /= (U**2).sum(axis=0)**0.5

    # sorted and binned neural activity
    X_embedding = bin1d(S[rm.isort], bin_size)
    rm.embedding_clust.max()


    # pupil resampled at neural times with some smoothing
    pupil = prediction.resample_frames(dat['behaviors']['pupil_moves'][np.newaxis,:], 
                            dat['timing']['cam_times'], 
                            dat['timing']['neural_times'], 
                           ).T
    
    S_pred, latents, itrain, itest, ival, state_dict = prediction.fit_neurons(dat, U[:,:min(U.shape[-1],128)], 
                                                                                device=device)

    
    S_pred = (S_pred * S_std) + S_mean
    torch.save(state_dict, root + f'proc/net_state_dict_{mouse}_D{day}')

    ve = 1 - ((S_pred[:,itest] - S[:,itest])**2).mean(axis=1) / ((S[:,itest] - S[:,itest].mean())**2).mean(axis=1)
    ve2 = 1 - ((S_pred[:,itest] - S[:,itest])**2).mean() / ((S[:,itest] - S[:,itest].mean())**2).mean()
    print('> TEST VAREXP: ', ve.mean(), ve2)
    ve = 1 - ((S_pred[:,ival] - S[:,ival])**2).mean(axis=1) / ((S[:,ival] - S[:,ival].mean())**2).mean(axis=1)
    ve2 = 1 - ((S_pred[:,ival] - S[:,ival])**2).mean() / ((S[:,ival] - S[:,ival].mean())**2).mean()
    print('> VAL VAREXP: ', ve.mean(), ve2)


    latents_itest = latents[itest].reshape(-1, latents.shape[-1])
    dists = ((latents_itest[:,:,np.newaxis] - latents[itrain].T)**2).sum(axis=1)
    imins = dists.argsort(axis=1)

    n_neighbors = 50
    S_pred_test_neigh = S[:,itrain][:,imins[:,:n_neighbors]].mean(axis=-1)

    ve = 1 - ((S_pred_test_neigh - S[:,itest])**2).mean(axis=1) / ((S[:,itest] - S[:,itest].mean())**2).mean(axis=1)
    ve2 = 1 - ((S_pred_test_neigh - S[:,itest])**2).mean() / ((S[:,itest] - S[:,itest].mean())**2).mean()
    print('> TEST NEIGHBOR VAREXP: ', ve.mean(), ve2)

    clu_id = dat['neurons']['clu_id']
    probe_id = dat['neurons']['probe_id']
    
    proc = {'mouse': mouse,  ## mouse name
            'day': day, ## recording day
            'isort': rm.isort, ## sorting by rastermap for visualization
            'sareas': sareas, ## number of area for each unit: run area_names[sareas] to get a name for each unit
            'area_names': area_names, ## area names corresponding to numbers in sarea
            'S_test': S[:,itest], ## spiking data at 10Hz during test periods
            'S_val': S[:, ival], ## spiking data at 10Hz during validation periods
            'S_pred_test': S_pred[:,itest], ## prediction from behavior of spiking data at 10Hz during test periods
            'S_pred_test_neigh': S_pred_test_neigh, ## prediction from nearest neighbors of behavior of spiking data at 10Hz during test periods
            'S_pred_val': S_pred[:,ival], ## prediction from behavior of spiking data at 10Hz during validation periods
            'latents': latents,
            'pupil': pupil, ## pupil abs change sampled at 10Hz
            'itest': itest, ## timepoints of test frames at 10Hz
            'ival': ival, ## timepoints of validation frames at 10Hz
            'clu_id': clu_id, ## cluster ids for the units from "clu"
            'probe_id': probe_id, ## probe ids for the units (0 or 1)
            'laser_durations': dat['timing']['laser_durations'],
            'laser_hz': dat['timing']['laser_hz']
            }

    np.save(root + f'proc/procnn_{mouse}_D{day}.npy', proc)
    
def get_stats(root):
    syncope = mat73.loadmat(root + 'Data/syncope_bouts_fix.mat')
    fnames = np.array([f[0] for f in syncope['matfiles']['name']])
    
    firing_rates = np.zeros((0,))
    firing_variances = np.zeros((0,))
    rasters_summary = np.zeros((0,120*10), 'int')
    latenciesS = np.zeros((0,))
    rastersS = np.zeros((0,2500), 'int')
    latenciesO = np.zeros((0,))
    rastersO = np.zeros((0,2500), 'int')
    latenciesC = np.zeros((0,))
    rastersC = np.zeros((0,2500), 'int')
    latenciesS0 = np.zeros((0,))
    rastersS0 = np.zeros((0,2500), 'int')
    latenciesR = np.zeros((0,))
    rastersR = np.zeros((0,2500), 'int')
    latenciesL = np.zeros((0,))
    rastersL = np.zeros((0,239), 'int')
    suppresseds = np.zeros((0,))
    areas = np.zeros((0,))
    timeoffs = np.zeros((0,))
    laser_frs = np.zeros((0, 4), 'float32')


    for i, mouse in enumerate(mice):
        for day in days[i]:
            out = summary.summary_stats(root, mouse, day, (syncope, fnames))
            firing_rate, firing_variance, raster_summary, suppressed, laser_fr, latencyS, rasterS, latencyO, rasterO, latencyC, rasterC, latencyS0, rasterS0, latencyR, rasterR, timeoff, latencyL, rasterL, brain_areas = out

            firing_rates = np.append(firing_rates, firing_rate, axis=0)
            firing_variances = np.append(firing_variances, firing_variance, axis=0)
            rasters_summary = np.append(rasters_summary, raster_summary, axis=0)
            suppresseds = np.append(suppresseds, suppressed, axis=0)
            latenciesS = np.append(latenciesS, latencyS, axis=0)
            rastersS = np.append(rastersS, rasterS, axis=0)
            latenciesO = np.append(latenciesO, latencyO, axis=0)
            rastersO = np.append(rastersO, rasterO, axis=0)
            latenciesC = np.append(latenciesC, latencyC, axis=0)
            rastersC = np.append(rastersC, rasterC, axis=0)
            latenciesS0 = np.append(latenciesS0, latencyS0, axis=0)
            rastersS0 = np.append(rastersS0, rasterS0, axis=0)
            latenciesR = np.append(latenciesR, latencyR, axis=0)
            rastersR = np.append(rastersR, rasterR, axis=0)
            timeoffs = np.append(timeoffs, timeoff, axis=0)
            latenciesL = np.append(latenciesL, latencyL, axis=0)
            rastersL = np.append(rastersL, rasterL, axis=0)
            areas = np.append(areas, brain_areas)
            laser_frs = np.append(laser_frs, laser_fr, axis=0)
        
    np.savez('summary_stats.npz', latenciesS, rastersS, latenciesO, rastersO, latenciesC, rastersC, latenciesL, rastersL,
                            suppresseds, laser_frs, areas, firing_rates, firing_variances, rasters_summary, timeoffs)
    

def get_residuals(root):
    
    syncope = mat73.loadmat(root + 'Data/syncope_bouts_fix.mat')
    fnames = np.array([f[0] for f in syncope['matfiles']['name']])
    
    firing_rates = np.zeros(0)
    laser_res = np.zeros(0)
    syncon_res = np.zeros(0)
    syncoff_res = np.zeros(0)
    pre_res = np.zeros(0)
    areas = np.zeros((0,))

    rasters_res = np.zeros((0, 3, 1200), 'float32')
    rasters_S = np.zeros((0, 3, 1200), 'float32')
    laser_res = np.zeros(0)
    syncon_res = np.zeros(0)
    syncoff_res = np.zeros(0)
    lasoff_res = np.zeros(0)
    pre_res = np.zeros(0)

    for i, mouse in enumerate(mice):
        for day in days[i]:

            dat = np.load(root + f'proc/dat_{mouse}_D{day}.npy', allow_pickle=True).item()
            proc = np.load(root + f'proc/procnn_{mouse}_D{day}.npy', allow_pickle=True).item()
            inds = np.nonzero(np.array([mouse in f and ('D'+str(day) in f or 'Day'+str(day) in f) for f in fnames]))[0]
            syncope_on = syncope['thresh'][:,inds[-1]]
            syncope_off = syncope['thresh2'][:,inds[-1]]

            S = dat['neurons']['S']
            n_neurons = S.shape[0]
            S_mean = S.mean(axis=1).copy()[:,np.newaxis]
            S_std = S.std(axis=1).copy()[:,np.newaxis]
            print('S = ', S.shape)

            firing_rate = (proc['S_val'] * 10).mean(axis=1)
            firing_variance = proc['S_val'].var(axis=1)

            splits = prediction.split_batches(dat['timing']['laser_times_n'], 
                                      dat['timing']['cam_times'], 
                                      dat['timing']['neural_times'], 
                                      l_train=500
                                     )
            itrain, ival, itest = splits[:3]
            n_neurons = S.shape[0]
            n_stims = len(itest) // 610
            if proc['laser_hz'][-1] == 20:
                inds = np.arange(0, n_stims)[-4:]
            else:
                inds = np.arange(0, n_stims)[1:-1]

            laser_hz = proc['laser_hz'][inds]
            raster_S = np.zeros((n_neurons, 3, 1200), 'float32')
            raster_res = np.zeros((n_neurons, 3, 1200), 'float32')

            X = dat['behaviors']['movSVD'][:,:250].copy()
            X -= X.mean(axis=0)
            X /= X[:,0].std() / 2
            X2 = dat['behaviors']['motSVD'][:,:250].copy()
            X2 -= X2.mean(axis=0)
            X2 /= X2[:,0].std()
            X = np.hstack((X, X2))
            device = torch.device('cuda')
            state_dict = torch.load(root + f'proc/net_state_dict_{mouse}_D{day}')
            n_out = state_dict['features.out.weight'].shape[0]
            net = prediction.Wavelets(n_in=X.shape[1], n_kp=25, kernel_size=51, n_latent=100,
                                          n_out=n_out).to(device)
            net.load_state_dict(state_dict)

            for ii in range(3):

                i0 = itest.reshape(-1,610)[inds[ii],0] + 9
                tinds = np.arange(i0 - 30 * 10, i0 + 10 * 90)
                print(laser_hz[ii], tinds[0], tinds[-1])
                Si = S[:, tinds]
                raster_S[:, ii] = Si

                cam_times = dat['timing']['cam_times'] 
                neural_times = dat['timing']['neural_times']
                icam0 = np.abs(cam_times - neural_times[tinds[0]]).argmin()
                icam1 = np.abs(cam_times - neural_times[tinds[-1]]).argmin()
                inds_cam = np.arange(icam0, icam1+1)
                i_sample = np.abs(cam_times - neural_times[tinds][:,np.newaxis]).argmin(axis=-1)
                i_sample -= icam0

                Xi = torch.from_numpy(X[inds_cam]).unsqueeze(0)
                with torch.no_grad():
                    latents_window = net(Xi.to(device), i_sample)[1].cpu().numpy().squeeze()
                y_pred, latents = prediction.compute_preds_latents(net, X, splits)

                dists = ((latents_window[:,:,np.newaxis] - latents[itrain].T)**2).sum(axis=1)
                imins = dists.argsort(axis=1)

                n_neighbors = 50
                S_pred_neigh = S[:,itrain][:,imins[:,:n_neighbors]].mean(axis=-1)

                raster_res[:,ii] = Si.copy() - S_pred_neigh

                if ii==2:
                    residual = raster_res[:,ii].copy()
                    l0 = 10 * 30
                    dtlas = 20
                    laser_on = residual[:, l0 : l0+dtlas].mean(axis=-1) #/ S[:,ival].mean(axis=-1)
                    dtlas = 20
                    laser_pre = residual[:, l0-2-dtlas : l0-2].mean(axis=-1)

                    dtlas = 25
                    s0 = l0 + int(np.round(syncope_on[0]*10)) + 5
                    sync_on = residual[:, s0 : s0+dtlas].mean(axis=-1)

                    dtlas = 25
                    s0 = l0 + int(np.round(syncope_off[0]*10)) + 5
                    sync_off = residual[:, s0 : s0+dtlas].mean(axis=-1)

                    #s0 = l0 + int(np.round(max(30.0, syncope_off.max()) * 10))
                    s0 = l0 + int(np.round(30.0*10)) - 25
                    dtlas = 50
                    laser_off = residual[:, s0 : s0+dtlas].mean(axis=-1)

                    print(s0+dtlas)

            laser_res = np.append(laser_res, laser_on, axis=0)
            lasoff_res = np.append(lasoff_res, laser_off, axis=0)
            pre_res = np.append(pre_res, laser_pre, axis=0)
            syncon_res = np.append(syncon_res, sync_on, axis=0)
            syncoff_res = np.append(syncoff_res, sync_off, axis=0)
            rasters_res = np.append(rasters_res, raster_res, axis=0)
            rasters_S = np.append(rasters_S, raster_S, axis=0)

            brain_areas = proc['area_names'][proc['sareas']]
            areas = np.append(areas, brain_areas)

    return laser_res, lasoff_res, pre_res, syncon_res, syncoff_res, rasters_res, rasters_S

def get_behaviors(root):
    pup_good = [[True, True], [False], [False, True], [True, True], [True, True], [False, False], [True, True]]

    whisk_pre = np.nan*np.zeros((13,))
    whisk_post = np.nan*np.zeros((13,))
    whisk_mean = np.nan*np.zeros((13,))
    pup_pre = np.nan*np.zeros((13,))
    pup_post = np.nan*np.zeros((13,))
    pup_mean = np.nan*np.zeros((13,))
    whisk_syncope = np.nan*np.zeros((13,3,360))
    whisk_laser = np.nan*np.zeros((13,3,360))
    cam_laser_time = np.zeros((13,3))
    syncope_time = np.zeros((13,))

    kk = 0
    for i, mouse in enumerate(mice):
        for day in days[i]:

            video_file = glob(root + f'videos/{mouse}*_D{day}*down2_proc_all.npy')
            video_file = video_file[0] if len(video_file)>0 else glob(root + f'videos/{mouse[:-1]}*_D{day}*_down2_proc_all.npy')[0]
            beh = np.load(video_file, allow_pickle=True).item()
            if 1:
                pcs = beh['motSVD'][1][:,:10].copy()
                signs = np.sign(skew(pcs, axis=0))
                pcs *= signs
                sig_baseline = 10
                win = 30 * 30
                Flow = gaussian_filter1d(pcs, sig_baseline, axis=0)
                Flow = minimum_filter1d(Flow,    win, axis=0)
                Flow = maximum_filter1d(Flow,    win, axis=0)
                whisk = ((pcs - Flow)**2).sum(axis=-1)**0.5
            else:
                whisk = (beh['motion'][1])
                avgmotion = beh['avgmotion'][0].reshape(beh['Ly'][0], beh['Lx'][0])
                avgwhisk = avgmotion[np.ix_(beh['rois'][1]['yrange'], beh['rois'][1]['xrange'])]
                print(whisk.min(), avgwhisk.sum())
                whisk += avgwhisk.sum()
            chin = (beh['motSVD'][2][:,:10]**2).sum(axis=-1)**0.5
            blink = beh['blink'][0]
            n_frames = len(blink)

            pup_x, pup_y, pup_area = np.nan * np.zeros((3, n_frames))
            if pup_good[i][day-1]:
                if 'pupil' in beh and len(beh['pupil']) > 0:
                    pup_area = beh['pupil'][0]['area'] #* 500/40)**2
                    pup_y, pup_x = beh['pupil'][0]['com'].T
                    pup_x = pup_x[:,np.newaxis]
                    pup_y = pup_y[:,np.newaxis]
                    iblink = np.isnan(pup_area)
                else:

                    pupil_file = glob(root + f'Data/Pupil/*{mouse}_*D*{day}_*.mat')
                    from scipy.io import loadmat
                    try:
                        pupil = mat73.loadmat(pupil_file[0])
                    except:
                        pupil = loadmat(pupil_file[0])

                    pup_x = pupil['pupildata'][9:, 1::3]
                    pup_y = pupil['pupildata'][9:, 2::3]
                    pup_llh = pupil['pupildata'][9:, 3::3]
                    n_kp = pup_x.shape[-1]

                    filt_x = np.array([median_filter(pup_x[:,i], 5) for i in range(n_kp)]).T
                    filt_y = np.array([median_filter(pup_y[:,i], 5) for i in range(n_kp)]).T

                    llh_cutoff = 0.75
                    pup_area = np.pi * np.prod(((pup_x[:,::2] - pup_x[:,1::2])**2 + (pup_y[:,::2] - pup_y[:,1::2])**2)**0.5, axis=1)
                    pup_area /= 4 # used diameter, need radius

                    iblink = ((pup_llh < llh_cutoff)==0).sum(axis=1) == 0


                iblink_starts = np.nonzero(blink < 3000)[0]
                iblinks = np.unique((iblink_starts[:,np.newaxis] + np.arange(0, 5)).flatten())
                iblink[iblinks] = True

                pup_x[iblink] = np.nan
                pup_y[iblink] = np.nan
                pup_area[iblink] = np.nan
                
                
            try:
                dat = np.load(root + f'proc/dat_{mouse}_D{day}.npy', allow_pickle=True).item()
                laser_hz = dat['timing']['laser_hz']
                laser_times_n = dat['timing']['laser_times_n']
            except:
                from scipy.io import loadmat
                digArray = loadmat(root + 'Data/NPY2R_Ai32_18N_Day2_g0_imec1_digArray.mat')['digArray']
                # times of TTL pulses / rate of signal
                cam_times = np.nonzero(np.diff(digArray[1,:])==1)[0] / 12500
                # subsetting due to dropped frames
                cam_times = cam_times[np.linspace(0, len(cam_times)-1, n_frames).astype(int)]
                laser_times = np.unique(np.nonzero((np.diff(digArray[2,:])>0))[0])
                laser_times = np.hstack((laser_times[0], laser_times[1:][np.diff(laser_times) > 5 * 12500/10])).astype('float64')
                laser_times /= 12500
                laser_hz = [5, 10, 20, 20]
                laser_times_n = None
            
            inds = np.nonzero(np.array([mouse in f and ('D'+str(day) in f or 'Day'+str(day) in f) for f in fnames]))[0]
            syncope_on = syncope['thresh'][:,inds[-1]]
            syncope_off = syncope['thresh2'][:,inds[-1]]
            
            n_stims = len(laser_hz)
            if laser_hz[-1] == 20:
                inds = np.arange(0, n_stims)[-4:]
            else:
                inds = np.arange(0, n_stims)[1:-1]
            laser_hz = np.array(laser_hz)[inds]
            
            import torch
            splits = prediction.split_batches(dat['timing']['laser_times_n'], 
                                    dat['timing']['cam_times'], 
                                    dat['timing']['neural_times'], 
                                    l_train=500
                                    )
            itrain, ival, itest, itrain_cam = splits[:4]
            itrain_cam = itrain_cam[itrain_cam < min(n_frames, len(pup_area))]

            for istim in range(3):

                if laser_times_n is not None:
                    i0 = laser_times_n[inds[istim]] - 1
                    tinds = np.arange(i0 - 10 * 30, i0 + 10 * 33)
                    cam_times = dat['timing']['cam_times'] 
                    neural_times = dat['timing']['neural_times']
                    icam0 = np.abs(cam_times - neural_times[tinds[0]]).argmin()
                    icam1 = np.abs(cam_times - neural_times[tinds[-1]]).argmin()
                    icam_laser = np.abs(cam_times - neural_times[tinds[10*30]]).argmin()
                else:
                    i0 = laser_times[inds[istim]]
                    icam0 = np.abs(cam_times - (i0-30)).argmin()
                    icam1 = np.abs(cam_times - (i0+33)).argmin()
                    icam_laser = np.abs(cam_times - (i0)).argmin()
                inds_cam = np.arange(icam0, icam1+1)
            
                if istim==2:
                    son = syncope_on[0]
                else:
                    son = 10
                isyncope = icam_laser + int(np.round(30 * son))
                whisk_syncope[kk,istim] = whisk[isyncope - 180 : isyncope + 180]
                whisk_laser[kk,istim] = whisk[icam_laser - 180 : icam_laser + 180]
                cam_laser_time[kk,istim] = icam_laser

            sec = 5
            
            whisk_pre[kk] = np.nanmean(whisk[icam_laser - 30*sec : icam_laser])
            whisk_post[kk] = np.nanmean(whisk[icam_laser : icam_laser + 30*sec])
            whisk_mean[kk] = np.nanmean(whisk[itrain_cam])
            pup_pre[kk] = np.nanmean(pup_area[icam_laser - 30*sec : icam_laser])
            pup_post[kk] = np.nanmean(pup_area[icam_laser : icam_laser + 30*sec])
            pup_mean[kk] = np.nanmean(pup_area[itrain_cam])
            syncope_time[kk] = syncope_on[0] 
            
    return (whisk_pre, whisk_post, whisk_mean, pup_pre, pup_post, pup_mean, 
            whisk_syncope, whisk_laser, cam_laser_time, syncope_time)