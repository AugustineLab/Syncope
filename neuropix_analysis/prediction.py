import torch
from torch import nn
from torch.nn import functional as F
import numpy as np
import sys, time
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from collections import OrderedDict

def resample_frames(data, torig, tout):
    ''' resample data at times torig at times tout '''
    ''' data is components x time '''
    fs = torig.size / tout.size # relative sampling rate
    data = gaussian_filter1d(data, np.ceil(fs/4), axis=1)
    f = interp1d(torig, data, kind='linear', axis=-1, fill_value='extrapolate')
    dout = f(tout)
    return dout


def copy_state(model):
    """
    Given PyTorch module `model`, makes a copy of the state onto CPU.
    Args:
        model: PyTorch module to copy state dict of

    Returns:
        A copy of state dict with all tensors allocated on the CPU
    """
    copy_dict = OrderedDict()
    state_dict = model.state_dict()
    for k, v in state_dict.items():
        copy_dict[k] = v.cpu() if v.is_cuda else v.clone()
    return copy_dict
   
def gabor_wavelet(sigma, f, ph, n_pts=201, is_torch=False):
    x = np.linspace(0, 2*np.pi, n_pts+1)[:-1].astype('float32')
    cos = np.cos
    sin = np.sin
    exp = np.exp
    xc = x - x.mean()
    cosine = cos(ph + f * xc)
    gaussian = exp(-(xc**2) / (2*sigma**2))
    G = gaussian * cosine
    G /= (G**2).sum()**0.5
    return G


def split_batches(laser_times_n, tcam, tneural, l_train=500, l_test=610, itrain=None, itest=None):
    n_t = len(tneural)
    if itrain is None or itest is None:
        
        inds_test = laser_times_n - 10
        itest = np.zeros(len(tneural), 'bool')
        itest[(inds_test[:,np.newaxis] + np.arange(0, l_test, 1, int)).flatten()] = True
        
        inds_train = []
        i = 0
        while i < n_t - l_train:
            if itest[i : i + l_train].sum() == 0:
                inds_train.append(i)
            i += l_train
        inds_train = np.array(inds_train)
        inds_val = inds_train[::5]
        inds_train = inds_train[~np.isin(inds_train, inds_val)]
        
        itrain = (inds_train[:,np.newaxis] + np.arange(0, l_train, 1, int)).flatten()
        ival = (inds_val[:,np.newaxis] + np.arange(0, l_train, 1, int)).flatten()
        itest = np.nonzero(itest)[0]
        
    # find itrain and itest in cam inds
    f = interp1d(tcam, np.arange(0, len(tcam)), kind='nearest', axis=-1,
                fill_value='extrapolate', bounds_error=False)

    inds_cam_train = f(tneural[inds_train]).astype('int')
    inds_cam_val = f(tneural[inds_val]).astype('int')
    inds_cam_test = f(tneural[inds_test]).astype('int')

    l_cam_train = int(np.ceil(np.diff(tneural).mean() / np.diff(tcam).mean() * l_train))
    l_cam_test = int(np.ceil(np.diff(tneural).mean() / np.diff(tcam).mean() * l_test))

    # create itrain and itest in cam inds
    itrain_cam = (inds_cam_train[:,np.newaxis] + np.arange(0, l_cam_train, 1, int)).flatten()
    ival_cam = (inds_cam_val[:,np.newaxis] + np.arange(0, l_cam_train, 1, int)).flatten()
    itest_cam = (inds_cam_test[:,np.newaxis] + np.arange(0, l_cam_test, 1, int)).flatten()
    
    itrain_cam = np.minimum(len(tcam)-1, itrain_cam)
    ival_cam = np.minimum(len(tcam)-1, ival_cam)
    itest_cam = np.minimum(len(tcam)-1, itest_cam)

    # inds for downsampling itrain_cam and itest_cam
    itrain_sample = f(tneural[itrain]).astype(int)
    ival_sample = f(tneural[ival]).astype(int)
    itest_sample = f(tneural[itest]).astype(int)
    
    #print('checking samples in cam inds')
    #print(np.isin(itrain_sample, itrain_cam).sum(), len(itrain_sample))
    #print(np.isin(itest_sample, itest_cam).sum(), len(itest_sample))

    it = np.zeros(len(tcam), 'bool')
    it[itrain_sample] = True
    itrain_cam_unique = np.unique(itrain_cam)
    print(len(itrain_cam), len(itrain_cam_unique))
    itrain_sample = it[itrain_cam_unique].nonzero()[0]

    it = np.zeros(len(tcam), 'bool')
    it[ival_sample] = True
    ival_cam_unique = np.unique(ival_cam)
    print(len(ival_cam), len(ival_cam_unique))
    ival_sample = it[ival_cam].nonzero()[0]
    
    it = np.zeros(len(tcam), 'bool')
    it[itest_sample] = True
    itest_sample = it[itest_cam].nonzero()[0]

    print(len(itrain_sample), len(itest_sample))
    return itrain, ival, itest, itrain_cam, ival_cam, itest_cam, itrain_sample, ival_sample, itest_sample, l_cam_train, l_cam_test


class Wavelets(nn.Module):
    def __init__(self, n_in=28, n_kp=30, n_filt=10, kernel_size=201, n_out=128, n_latent=100):
        super().__init__()
        self.n_in = n_in
        self.n_kp = n_kp
        self.n_filt = n_filt
        self.n_out = n_out
        self.n_latent = n_latent
        self.kernel_size = kernel_size

        self.features = nn.Sequential()

        # combine keypoints into n_kp features
        self.features.add_module('linear0_0', nn.Linear(n_in//2, n_kp))
        self.features[-1].weight.data = 0.0001 * torch.randn(self.features[-1].weight.data.shape)
        
        self.features.add_module('linear0_1', nn.Linear(n_in//2, n_kp))
        self.features[-1].weight.data = 0.0001 * torch.randn(self.features[-1].weight.data.shape)
        
        
        # compute n_filt wavelet features of each one => n_filt * n_kp features
        self.features.add_module('wavelet0_0', nn.Conv1d(1, n_filt, kernel_size=kernel_size,
                                                      padding=kernel_size//2, bias=False))
        # initialize filters with gabors
        f = np.geomspace(1, 10, n_filt//2).astype('float32')
        gw0 = gabor_wavelet(10, f[:,np.newaxis], np.pi, n_pts=kernel_size)
        gw1 = gabor_wavelet(10, f[:,np.newaxis], np.pi/2, n_pts=kernel_size)
        #self.features[-1].weight.data = torch.from_numpy(gw0).unsqueeze(1)
        self.features[-1].weight.data = torch.from_numpy(np.vstack((gw0, gw1))).unsqueeze(1)
        
        # compute n_filt wavelet features of each one => n_filt * n_kp features
        self.features.add_module('wavelet0_1', nn.Conv1d(1, n_filt, kernel_size=kernel_size,
                                                      padding=kernel_size//2, bias=False))
        self.features[-1].weight.data = torch.from_numpy(np.vstack((gw0, gw1))).unsqueeze(1)
        
        # latent linear layer
        self.features.add_module('linear1', nn.Linear(n_filt * n_kp * 2 , n_latent))
        self.features[-1].weight.data = 0.0001 * torch.randn(self.features[-1].weight.data.shape)
        
        # tried filtering again, didn't help
        #self.features.add_module('wavelet1', nn.Conv1d(100, 100, kernel_size=5,
        #                                               groups=100,
        #                                              padding=5//2, bias=False))
        
        # output linear layer
        self.features.add_module('out', nn.Linear(n_latent, n_out)) #n_filt * n_kp, n_out)) #
        self.features[-1].weight.data = 0.0001 * torch.randn(self.features[-1].weight.data.shape)
        
        self.bias = nn.Parameter(torch.zeros((n_out)))
                                                       
    def forward(self, x, sample_inds):
        """ x is (n_batches, time, features)
            sample_inds is (sub_time) over batches
        """
        n_batches = x.shape[0]
        n_f = x.shape[-1]
        n_kp = self.n_kp
        n_filt = self.n_filt
        # x is (n_batches, time, features)
        out = self.features[0](x[..., :n_f//2]).transpose(2,1)
        # out is now (n_batches, n_kp, time)
        out = out.reshape(-1, out.shape[-1]).unsqueeze(1)
        # out is now (n_batches * n_kp, 1, time)
        out = self.features[2](out)  #[:,:,1+5:-self.kernel_size+5]
        # out is now (n_batches * n_kp, n_filt, time)
        out0 = out.reshape(n_batches, n_kp * n_filt, -1)
        # out0 now (n_batches, n_kp * n_filt, time)
        
        out = self.features[1](x[..., n_f//2:]).transpose(2,1)
        # out is now (n_batches, n_kp, time)
        out = out.reshape(-1, out.shape[-1]).unsqueeze(1)
        # out is now (n_batches * n_kp, 1, time)
        out = self.features[3](out)  #[:,:,1+5:-self.kernel_size+5]
        # out is now (n_batches * n_kp, n_filt, time)
        out1 = out.reshape(n_batches, n_kp * n_filt, -1)
        # out1 now (n_batches, n_kp * n_filt, time)
        
        out = torch.cat((out0, out1), dim=1)
        # out is now (n_batches, n_kp * n_filt * 2, time)
        out = out.transpose(2,1)
        
        out = out.reshape(-1, self.n_kp * self.n_filt * 2)
        out = out[sample_inds]
        out = F.relu(out)
        
        # out => (n_batches, downsampled_time, n_kp * n_filt)
        
        # linear layer
        wav = self.features[-2](out.reshape(n_batches, -1, out.shape[-1]))
        #wav = F.relu(out)
        
        # final layer
        out = self.features[-1](wav)
        #out = F.relu(out) + self.bias
        
        return out, wav


def train_epoch(net, optimizer, X_train, Y_train, itrain_sample_b, epoch=0, batch_size=1, smoothing_penalty=1.0, device=torch.device('cuda')):
    """ train epoch of wavelets model
    
    always batch_size = 1 (one segment of time) 
    
    """
    net.train()
    n_batches = X_train.shape[0]
    np.random.seed(epoch)
    rperm = np.random.permutation(n_batches)
    train_loss = 0
    gmax = 0
    for n in rperm:
        y_pred, wav = net(X_train[n].unsqueeze(0).to(device), itrain_sample_b[n])
        loss = ((y_pred - Y_train[n].unsqueeze(0).to(device))**2).mean()
        loss += smoothing_penalty * (torch.diff(net.features[1].weight)**2).sum()
        
        optimizer.zero_grad()
        loss.backward()

        gnorm = nn.utils.clip_grad_norm_(net.parameters(), 10)
        #print(f'{gnorm.item():.3f}', end='\t')
        gmax = max(gmax, gnorm)
        
        optimizer.step()
        train_loss += loss.item()
    train_loss /= n_batches

    return train_loss

def train_wavelets(net, X, Y, splits, learning_rate=1e-2, weight_decay=5e-6, smoothing_penalty=1.0, n_epochs=300, device=torch.device('cuda')):
    ### making segments as batches so that filters don't run across
    itrain, ival, itest, itrain_cam, ival_cam, itest_cam, itrain_sample, ival_sample, itest_sample, l_cam_train, l_cam_test = splits
    
    X_train = torch.from_numpy(X[itrain_cam].reshape(-1, l_cam_train, X.shape[-1])).float()
    Y_train = torch.from_numpy(Y[itrain]).float()
    Y_train = Y_train.reshape(X_train.shape[0], -1, Y.shape[-1])

    X_val = torch.from_numpy(X[ival_cam].reshape(-1, l_cam_train, X.shape[-1])).float()
    Y_val = torch.from_numpy(Y[ival]).float()
    Y_val = Y_val.reshape(X_val.shape[0], -1, Y.shape[-1])
    
    X_test = torch.from_numpy(X[itest_cam].reshape(-1, l_cam_test, X.shape[-1])).float()
    Y_test = torch.from_numpy(Y[itest]).float()
    Y_test = Y_test.reshape(X_test.shape[0], -1, Y.shape[-1])

    itrain_sample_b = itrain_sample.reshape(len(itrain_cam) // l_cam_train, -1).copy()
    itrain_sample_b -= itrain_sample_b[:,:1]
    
    #itest_sample_b = itest_sample.reshape(len(itest_cam) // l_cam_test, -1)
    
    n_periods = 3
    for i_period in range(n_periods):
        lr = learning_rate / (3 ** (i_period))
        print(lr)

        restore = (i_period > 0)
        if restore:
            net.load_state_dict(best_state_dict)
        else:
            varexp_max = 0

        optimizer = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=weight_decay)
        tic = time.time()

        n_epochs = n_epochs if i_period == 0 else 30

        for epoch in range(n_epochs):
            train_loss = train_epoch(net, optimizer, X_train, Y_train, itrain_sample_b, epoch=epoch,
                                    smoothing_penalty=smoothing_penalty, device=device)
            
            with torch.no_grad():
                net.eval()
                y_pred_val = net(X_val.to(device), ival_sample)[0].cpu()
                tl = ((y_pred_val - Y_val)**2).mean()
                vev = 1 - tl / (Y_val**2).mean() 
                if vev > varexp_max:
                    best_state_dict = copy_state(net)
                    varexp_max = vev

            if epoch%10==0 or epoch==n_epochs-1:
                with torch.no_grad():
                    y_pred_test = net(X_test.to(device), itest_sample)[0].cpu()
                    tl = ((y_pred_test - Y_test)**2).mean()
                    vet = 1 - tl / (Y_test**2).mean()   
                    print(f'{epoch}, train loss {train_loss:.4f}, varexp val {vev.item():.4f}, varexp test {vet.item():.4f}, time {time.time()-tic:.2f}s')
    
    # return y_pred_val and y_pred_test from best network
    net.load_state_dict(best_state_dict)
    with torch.no_grad():
        net.eval()
        y_pred_val = net(X_val.to(device), ival_sample)[0].cpu()
        y_pred_test = net(X_test.to(device), itest_sample)[0].cpu()
        
    return y_pred_test.detach().cpu().numpy(), y_pred_val.detach().cpu().numpy(), best_state_dict


def compute_preds_latents(net, X, splits, device=torch.device('cuda')):
    """ get latents from net in response to X """
    itrain, ival, itest, itrain_cam, ival_cam, itest_cam, itrain_sample, ival_sample, itest_sample, l_cam_train, l_cam_test = splits
    
    X_train = torch.from_numpy(X[itrain_cam].reshape(-1, l_cam_train, X.shape[-1])).float()
    X_val = torch.from_numpy(X[ival_cam].reshape(-1, l_cam_train, X.shape[-1])).float()
    X_test = torch.from_numpy(X[itest_cam].reshape(-1, l_cam_test, X.shape[-1])).float()

    itrain_sample_b = itrain_sample.reshape(len(itrain_cam) // l_cam_train, -1).copy()
    
    n_time = max(itrain.max(), ival.max(), itest.max()) + 1

    with torch.no_grad():
        y_pred_itrain0, wav_itrain0 = net(X_train[:12].to(device), itrain_sample_b[:12].flatten())
        y_pred_itrain1, wav_itrain1 = net(X_train[12:].to(device), itrain_sample_b[12:].flatten() - itrain_sample_b[12,0])
        wav_itrain = np.vstack((wav_itrain0.cpu().numpy(), wav_itrain1.cpu().numpy()))
        y_pred_itrain = np.vstack((y_pred_itrain0.cpu().numpy(), y_pred_itrain1.cpu().numpy()))
        y_pred_ival, wav_ival = net(X_val.to(device), ival_sample)
        wav_ival = wav_ival.cpu().numpy()
        y_pred_itest, wav_itest = net(X_test.to(device), itest_sample)
        wav_itest = wav_itest.cpu().numpy()

    y_pred = np.zeros((n_time, y_pred_itrain.shape[-1]), 'float32')
    y_pred[itrain] = y_pred_itrain.reshape(-1, y_pred_itrain.shape[-1])
    y_pred[ival] = y_pred_ival.reshape(-1, y_pred_itrain.shape[-1]).cpu().numpy()
    y_pred[itest] = y_pred_itest.reshape(-1, y_pred_itrain.shape[-1]).cpu().numpy()
    latents = np.zeros((n_time, wav_itrain.shape[-1]), 'float32')
    latents[itrain] = wav_itrain.reshape(-1, wav_itrain.shape[-1])
    latents[ival] = wav_ival.reshape(-1, wav_itrain.shape[-1])
    latents[itest] = wav_itest.reshape(-1, wav_itrain.shape[-1])

    return y_pred, latents    
                

def fit_neural_pcs(dat, U, state_dict=None, device=torch.device('cuda')): 
    """ predict PCs of neural activity using movSVD and motSVD """
           
    splits = split_batches(dat['timing']['laser_times_n'], 
                          dat['timing']['cam_times'], 
                          dat['timing']['neural_times'], 
                          l_train=500
                         )
    itrain, ival, itest = splits[:3]

    Xs = []
    for key in ['movSVD', 'motSVD']:
        X = dat['behaviors'][key][:,:250].copy()
        X -= X.mean(axis=0)
        X /= X[:,0].std()
        if key=='movSVD':
            # (weight movSVD less than motSVD)
            X /= 2
        Xs.append(X)
    X = np.hstack(Xs)
    print(X.shape)
    
    Y = dat['neurons']['S'].T @ U

    Yt_mean = 0 * Y.mean(axis=1)
    Y -= Yt_mean[:,np.newaxis]
    Y_mean = Y.mean(axis=0)
    Y -= Y_mean
    
    net = Wavelets(n_in=X.shape[1], n_kp=25, kernel_size=51, n_latent=100,
                              n_out=Y.shape[1]).to(device)
    net.features[2].weight.requires_grad = False
    net.features[3].weight.requires_grad = False
    if state_dict is None:
        y_pred_test, y_pred_val, state_dict = train_wavelets(net, X, Y, splits, 
                                                            n_epochs=200,
                                                            smoothing_penalty=1.0,
                                                            learning_rate=5e-4, 
                                                            weight_decay=5e-7, 
                                                            device=device)
            
    y_pred, latents = compute_preds_latents(net, X, splits)

    return y_pred, latents, Yt_mean, Y_mean, itrain, itest, ival, state_dict

def fit_neurons(dat, U, state_dict=None, device=torch.device('cuda')):
    v_pred, latents, vt_mean, v_mean, itrain, itest, ival, state_dict = fit_neural_pcs(dat, U, state_dict=state_dict, device=device)
    v_pred += v_mean
    
    S_pred = U @ v_pred.T
    
    return S_pred, latents, itrain, itest, ival, state_dict

def plot_prediction(fig, X_embedding, X_embedding_pred, inds, laser_times, pupil, vmin=0, vmax=1):
    
    ax = fig.add_subplot(3,1,1)
    ax.imshow(X_embedding[:,inds], aspect='auto', vmin=0, vmax=vmax, cmap='gray_r')
    n_feats = X_embedding.shape[0]
    for lt in laser_times:
        ax.plot([lt, lt], [0, n_feats], 'y')
    ax.set_ylim([0, n_feats])
    ax.set_title('neural activity')
    ax.set_ylabel('binned neurons')
    ax.set_xlabel('time @ 10Hz')
    ax.plot(pupil[inds] / pupil.max() * 100)
    
    ax = fig.add_subplot(3,1,2)
    ax.imshow(X_embedding_pred[:,inds], aspect='auto', cmap='gray_r', vmin=vmin, vmax=vmax)
    for lt in laser_times:
        ax.plot([lt, lt], [0, n_feats], 'y')
    ax.set_ylim([0, n_feats])
    ax.set_title('prediction from behavior')
    ax.set_ylabel('binned neurons')
    ax.set_xlabel('time @ 10Hz')
    ax.plot(pupil[inds] / pupil.max() * 100)
    
    ax = fig.add_subplot(3,1,3)
    ax.imshow(X_embedding[:,inds] - X_embedding_pred[:,inds], aspect='auto', cmap='bwr', vmin=-2, vmax=2)
    for lt in laser_times:
        ax.plot([lt, lt], [0, n_feats], 'y')
    ax.set_ylim([0, n_feats])
    ax.set_title('residual')
    ax.set_ylabel('binned neurons')
    ax.set_xlabel('time @ 10Hz')
    ax.plot(pupil[inds] / pupil.max() * 100)
    
    
    














