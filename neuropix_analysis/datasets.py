
from glob import glob
import mat73
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import zscore
from scipy.sparse import csr_matrix
from scipy.ndimage import gaussian_filter1d
from natsort import natsorted
from rastermap.utils import bin1d

def load_neurons_behaviors(root, mouse, day):

    files = natsorted(glob(root + f'Data/*{mouse}*_D{day}*_imec*'))
    files = files if len(files) > 0 else natsorted(glob(root + f'Data/*{mouse}*_Day{day}_g0_imec*'))
    print(files)
    dats = [mat73.loadmat(f) for f in files]
    st = np.hstack(tuple(d['sp']['st'] for d in dats))
    nc0 = int(dats[0]['sp']['clu'].max()+1)
    nc1 = 0 if len(dats)==1 else int(dats[1]['sp']['clu'].max()+1)
    print('units per probe = ', nc0, nc1)
    clu = np.hstack(tuple(d['sp']['clu'] if k==0 else d['sp']['clu'] + nc0
                          for k,d in enumerate(dats)))
    probe_id = np.hstack((np.zeros((nc0,),'int'), np.ones((nc1,), 'int')))
    clu_id = np.hstack((np.arange(0, nc0, 1, 'int'), np.arange(0, nc1, 1, 'int')))
    brain_areas = np.hstack(tuple(np.array(d['cluster_region'], dtype='object')[2] for d in dats))
    digArray = dats[0]['digArray']
    TrainInfo = dats[0]['TrainInfo']

    video_file = glob(root + f'videos/{mouse}*_D{day}*down2_proc_all.npy')
    video_file = video_file[0] if len(video_file)>0 else glob(root + f'videos/{mouse[:-1]}*_D{day}*_down2_proc_all.npy')[0]
    proc = np.load(video_file, allow_pickle=True).item()
    motSVD = proc['motSVD'][1]
    movSVD = proc['movSVD'][1]
    n_frames = motSVD.shape[0]

    pupil_file = glob(root + f'Data/Pupil/{mouse}_*D{day}*_pupil.mat')
    if len(pupil_file) > 0:
        pupil = mat73.loadmat(pupil_file[0])
        pupil_moves = pupil['pupildata'][:,1:]
        pupil_moves -= pupil_moves.mean(axis=0)
        pupil_moves = np.abs(pupil_moves).sum(axis=-1)
    else:
        pupil_moves = np.zeros(n_frames, 'float32')

    # times of TTL pulses / rate of signal
    cam_times = np.nonzero(np.diff(digArray[1,:])==1)[0] / 12500
    # subsetting due to dropped frames
    cam_times = cam_times[np.linspace(0, len(cam_times)-1, n_frames).astype(int)]
    pupil_moves = pupil_moves[np.linspace(0, len(pupil_moves)-1, n_frames).astype(int)]

    print(len(cam_times), motSVD.shape)

    # bin in 100Hz (assuming only 1 spike per 10ms)
    srate = 100

    S = csr_matrix((np.ones(len(st), 'uint8'), 
                   (clu, np.round(st * srate).astype(int))), 
                   shape=(clu.max()+1, 
                          np.round(st * srate).astype(int).max()+1))

    # discard cells < 0.25 Hz
    fr = S.mean(axis=1) * srate
    goodcells = np.array(fr > 0.25).flatten()
    if mouse=='59L' and day==2:
        fr0 = np.array(fr > 0.).flatten()
        S = S[fr0]
        goodcells = goodcells[fr0]
        clu_id = clu_id[fr0]
        probe_id = probe_id[fr0]
        #brain_areas = brain_areas[fr0]

    S = S[goodcells]
    areas = brain_areas[goodcells]
    clu_id = clu_id[goodcells]
    probe_id = probe_id[goodcells]
    areas = np.array([a[0] for a in areas], dtype='object')

    # bin in 10Hz
    binsize = 10
    srate /= binsize
    ns = (S.shape[1]//binsize)*binsize
    nn = S.shape[0]
    nt = ns//binsize
    Sb = np.zeros((nn, nt), 'uint8')
    for n in range(nn):
        Sb[n] = np.array(S[n,:ns].reshape(nt, binsize).sum(axis=-1)).flatten()

    # neural times (crop by camera times)
    neural_times = np.arange(0, nt/srate, 1/srate)
    tinit = np.nonzero(neural_times > cam_times[0])[0][0] + 1
    tend = np.nonzero(neural_times < cam_times[-1])[0][-1] - 1
    Sb = Sb[:, tinit : tend]
    neural_times = neural_times[tinit : tend]

    #laser_times = TrainInfo[2]
    laser_times = np.unique(np.nonzero((np.diff(digArray[2,:])>0))[0])
    laser_times = np.hstack((laser_times[0], laser_times[1:][np.diff(laser_times) > 5 * 12500/10]))
    laser_durations = np.round(TrainInfo[-1]).astype(int)
    laser_hz = np.round(TrainInfo[1]).astype(int)
    laser_times_n = np.floor(laser_times / 12500 * srate).astype(int)
    laser_times_n -= tinit


    # exclude laser times when detecting drifting neurons
    inds_test = laser_times_n - 10
    l_test = 900
    itest = np.zeros(len(neural_times), 'bool')
    itest[(inds_test[:,np.newaxis] + np.arange(0, l_test, 1, int)).flatten()] = True
    itrain = ~itest
    rates = bin1d(Sb[:, itrain].T, int(srate)).T
    smoothed = gaussian_filter1d(rates, 300, axis=1)
    drift_ratio = smoothed.max(axis=1) / (1e-1 + smoothed.min(axis=1))

    drift_threshold = 5.
    print((drift_ratio < drift_threshold).mean())
    S = Sb[drift_ratio < drift_threshold].copy().astype('float32')
    #S = zscore(S, axis=1)
    area_names, sareas = np.unique(areas[drift_ratio < drift_threshold], return_inverse=True)
    clu_id = clu_id[drift_ratio < drift_threshold]
    probe_id = probe_id[drift_ratio < drift_threshold]

    dat = {
           'neurons': {'S': S,
                       'sareas': sareas,
                       'area_names': area_names,
                       'probe_id': probe_id,
                       'clu_id': clu_id
                      },
           'behaviors': {'movSVD': movSVD,
                         'motSVD': motSVD,
                         'pupil_moves': pupil_moves,
                        },
           'timing': {'neural_times': neural_times,
                      'cam_times': cam_times,
                      'laser_times_n': laser_times_n,
                      'laser_durations': laser_durations,
                      'laser_hz': laser_hz
                     }
          }
    np.save(root + f'proc/dat_{mouse}_D{day}.npy', dat)
    
    return dat


def region_categories():
    mainregions = ['ACA','AI','BLA','BMA','BST','DORpm','DORsm','HPC','ILA','LA','LS',
                    'LZ','MBmot','MBsta','MEZ','MOp','MOs','OLF','ORB',
                    'p-mot','PA','PALd','PALm','PALv','PIR','PL','PVR','PVZ',
                    'RSP','SC','SSp','SSs','STRd','STRv','sAMY','VISC','PTLp','VS','Fiber','error']
    region=[]
    region.append(['ACAd1','ACAd2/3', 'ACAd5', 'ADAd6a','ACAv1','ACAv2/3','ACAv5','ACAv6a','ACAv6b', 'ACAd6a'])
    region.append(['AIp2/3', 'AIp5', 'AId6a','AIp6a', 'AIp6b', 'AIv6a','CLA','GU6a']) # added AIp2/3,5 ***
    region.append(['BLAa','BLAp', 'CTXsp'])
    region.append(['BMAa','BMAp'])
    region.append(['BST'])
    region.append(['AD','AMd','AMv','AV','CL','CM','Eth','IAM','IMD','LD','LH','LP','MD','MH','PCN','PF','PO','PT','PVT','RE','RH','RT','SMT','SubG','TH','Xi','IAD','SGN'])
    region.append(['LGd-co','LGd-ip','LGd-sh','LGv','VAL','VM','VPL','VPLpc','VPM','VPMpc','PoT','SPFp'])#DORsm (add PoT,SPFp)
    region.append(['CA1','CA2','CA3','DG-mo','DG-po','DG-sg','HPF','ProS','SUB','IG']) #HPC (add IG)
    region.append(['ILA1','ILA2/3','ILA5','ILA6a','ILA6b'])#}; %ILA
    region.append(['LA'])#}; %LA
    region.append(['LSc','LSr','LSv','SF'])#}; %LS
    region.append(['FF','HY','LHA','LPO','PSTN','ZI'])#}; %LZ
    region.append(['APN','MB','MRN','MT','NOT','PAG','PPT','PR','RN','RR','SNr','VTA'])#}; %MBmot
    region.append(['PPN','SNc'])#}; %MBsta
    region.append(['MPN','PH','PMv','PMd','TMv','TU','VMH','AHN'])#}; %MEZ (add AHN)
    region.append(['MOp2/3','MOp5','MOp6a','MOp6b'])#}; %MOp
    region.append(['MOs1','MOs2/3','MOs5','MOp6a','MOp6b'])#}; %MOs
    region.append(['AON','COApl','COApm','DP','EPd','EPv','NLOT3','OLF','PAA','TTd','TR'])#}; %OLF (add TR)
    region.append(['FRP5','FRP6a','ORBl1','ORB12/3','ORBl5','ORBl6a','ORBm1','ORBm2/3','ORBm5','ORBvl6a'])#}; %ORB
    region.append(['PG','PRNc','PRNr','TRN'])#}; %p-mot
    region.append(['PA'])#}; %PA
    region.append(['GPe','GPi','PAL'])#}; %PALd
    region.append(['MS','NDB','TRS'])#}; %PALm
    region.append(['Sl'])#}; %PALv
    region.append(['PIR'])#}; %PIR
    region.append(['PL1','PL2/3','PL5','PL6a','PL6b'])#}; %PL
    region.append(['DMH','SBPV'])#}; %PVR
    region.append(['PVH','PVHd','PVi'])#}; %PVZ
    region.append(['RSPv2/3','RSPv5','RSPv6a'])#}; %RSP
    region.append(['SCdg','SCdw','SCig','SCiw','SCop'])#}; %SC
    region.append(['SSp-n4','SSp-n5','SSp-n6a','SSp-ul6a','SSP-ul6b','SSp-bfd6a','SSp-bfd6b','SSp-ll2/3','SSp-ll4','SSp-ll5','SSp-ll6a','SSp-ll6b','SSp-n1','SSp-n2/3','SSp-tr1','SSp-ul1','SSp-ul2/3','SSp-ul4','SSp-ul5','SSp-ul6b','SSp-un2/3','SSp-un4','SSp-un5','SSp-un6a','SSp-un6b', 'SSp-m5', 'SSp-m6a'])#}; %SSp (add SSp-bfd6a->)
    region.append(['SSs5','SSs6a','SSs6b'])#}; %SSs
    region.append(['CP','FS','STR'])#}; %STRd
    region.append(['ACB','OT'])#}; %STRv
    region.append(['AAA','CEac','CEAl','CEAm','IA','MEA'])#}; %sAMY
    region.append(['VISC6a','VISC6b'])#}; %VISC
    region.append(['VISa1','VISa2/3','VISa4','VISa5','VISa6a','VISa6b'])#}; %PTLp
    region.append(['SEZ','V3','VL'])#}; %VS- ventricular system
    region.append(['aco','act','alv','amc','bsc','ccb','ccg','ccs','cing','cpd','ec','ee','fa','fi','fiber tracts','em','ml','or','vhc','scwm'])#}; %Fiber tracts
    region.append(['outside_brain','root','o'])#}; % channel not in brain
    
    superregions = ['cortex', 'amygdala', 'hippocampus', 'striatum', 'pallidum', 'hypothalamus', 'midbrain+thalamus', 'fibertracts', 'error']
    super_region = []
    super_region.append(['ACA', 'AI', 'ILA', 'MOp', 'MOs', 'OLF', 'ORB', 'PIR', 'PL', 'RSP', 'SSp', 'SSs', 'VISC', 'PTLp'])
    super_region.append(['BLA', 'BMA', 'LA', 'PA'])
    #super_region.append([])
    super_region.append(['HPC'])
    super_region.append(['LS', 'STRd', 'STRv', 'sAMY'])
    super_region.append(['PALd', 'PALm', 'PALv', 'BST'])
    super_region.append(['LZ', 'MEZ', 'PVR', 'PVZ'])
    super_region.append(['MBmot', 'MBsta', 'SC', 'DORpm', 'DORsm'])
    #super_region.append(['p-mot'])
    super_region.append(['VS', 'Fiber'])
    super_region.append(['error'])




    return mainregions, region, superregions, super_region

