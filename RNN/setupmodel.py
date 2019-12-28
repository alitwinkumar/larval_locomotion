#this file is part of zarin_et_al_multilayer_2019
#Copyright (C) 2019 Ashok Litwin-Kumar
#see README for more information

import tensorflow as tf
import numpy as np
import h5py

def loadconns():
    f = h5py.File("data.h5","r")
    #after loading, rows are postsynaptic and columns are presynaptic
    Jpm = ((f["Jpm"][:]).T).astype(np.float32)
    Jpp = ((f["Jpp"][:]).T).astype(np.float32)

    pnames = f["p"][:]
    pnames_orig = f["p_orig"][:]
    mnames = f["m"][:]
    types = f["nt"][:]
    mnorder = f["mnorder"][:]
    f.close()

    return Jpm,Jpp,pnames,pnames_orig,mnames,types,mnorder

def genpulse(Tstop,dt,pstart,pend,Trise):
    T = int(Tstop/dt)

    p = np.zeros(T)

    istart = int(pstart/dt)
    iend = int(pend/dt)+1
    irise = int(Trise/dt)

    p[istart:iend] = 1.
    p[istart:(istart+irise)] = np.sin(np.pi*np.arange(irise)/(2.*irise))
    p[(iend-irise):iend] = np.flipud(np.sin(np.pi*np.arange(irise)/(2.*irise)))
    return p

def gentargets(mnorder,B,S,dt):
    M = mnorder.shape[1]

    Tstop = 6
    Tpulserise = 1. #pulse length
    Tpulse = 2.
    dtpulse = 0.25 #time between start of one pulse and the next
    dtpulse_end = 0.125 #time between start of one pulse and the next
    tstart = 1.
    segdelay = 1.
    T = int(Tstop/dt)

    Npulse = np.max(mnorder)

    pstarts = np.zeros([2,M])
    pends = np.zeros([2,M])

    for mi in range(M):
        for bi in range(B):
            segoffset = 0
            if (bi == 1) and (mi >= int(M/2)): #A2 MNs fire later during backward
                segoffset = segdelay
            elif (bi == 0) and (mi < int(M/2)): #A1 MNs fire later during forward
                segoffset = segdelay

            pstarts[bi,mi] = tstart + (mnorder[bi,mi]-1)*dtpulse + segoffset
            pends[bi,mi] = tstart + Tpulse + (mnorder[bi,mi]-1)*dtpulse_end + segoffset

    mtarg = np.zeros([T,B,M])
    costw = np.zeros([T,B,M])
    s = np.zeros([T,B,S],dtype=np.float32)

    for mi in range(M):
        for bi in range(B):
            if mnorder[bi,mi] > 0:
                mtarg[:,bi,mi] = genpulse(Tstop,dt,pstarts[bi,mi],pends[bi,mi],(pends[bi,mi]-pstarts[bi,mi])/2.)
                Norder = np.sum(mnorder[bi,:] == mnorder[bi,mi])/2
                costw[:,bi,mi] = 1. / np.sqrt(Norder)
            else: #unspecified target
                costw[:,bi,mi] = 0.

    seg1inds = np.where(np.sum(mtarg[:,1,0:int(M/2)],1)>0)[0]
    seg2inds = np.where(np.sum(mtarg[:,1,int(M/2):M],1)>0)[0]

    square_pulse = np.zeros(T,dtype=np.float32)
    istart = int(tstart/dt)
    iend = np.max(seg2inds)
    square_pulse[istart:] = 1.

    for bi in range(int(B/2)):
        s[:,bi,0] = square_pulse

    for bi in range(int(B/2),B):
        s[:,bi,1] = square_pulse

    return mtarg,costw,s,seg1inds,seg2inds

def rectify(x):
    return x*(x>0)

def initJpp(Jpp0,types):
    J = np.zeros(Jpp0.shape,dtype=np.float32)
    N = J.shape[0]
    N2 = J.shape[1]
    J = np.copy(Jpp0)

    for qi in range(N):
        if types[qi] == 'inh':
            J[qi,:] = -J[qi,:]
        elif types[qi] == '':
            J[qi,:] = -J[qi,:]
    return J

Jpm0,Jpp0,pnames,pnames_orig,mnames,types,mnorder = loadconns()

Nepochs = 1000
alpha = np.power(np.arange(Nepochs)/(Nepochs),2)
alpha_J = .1
alpha_seg = .1
alpha1827 = 0.05

M = len(mnames) #number of MNs
N = len(pnames) #number of PMNs
S = 2
B = 2
Ncycles = 1
dt = 0.05
lr = np.logspace(-2,-3,Nepochs).astype(np.float32)
taumin = .05
taumax = 1.
tau0 = 0.2
g0 = 1.
wsp0 = 0.2 + 0.1*np.random.rand(S,N) #strength of pulse at initial state
sigp0 = .2 #variability of initial state

ind18 = np.where(['A18b_' in s for s in pnames])[0]
ind27 = np.where(['A27h_' in s for s in pnames])[0]

mtarg,costw,s,seg1inds,seg2inds = gentargets(mnorder,B,S,dt)

T = mtarg.shape[0]

def scan_step(prevstate,inputs):
    with tf.variable_scope('model',reuse=True):
        Jpp = tf.get_variable('Jpp',shape=[N,N])
        Jpp_mask = tf.get_variable('Jpp_mask',shape=[N,N])
        Jpm = tf.get_variable('Jpm',shape=[N,M])
        Jpm_mask = tf.get_variable('Jpm_mask',shape=[N,M])
        bm = tf.get_variable('bm',shape=[M])
        bp = tf.get_variable('bp',shape=[N])
        taum = tf.get_variable('taum',shape=[M])
        taup = tf.get_variable('taup',shape=[N])
        gm = tf.get_variable('gm',shape=[M])
        gp = tf.get_variable('gp',shape=[N])
        wsp = tf.get_variable('wsp',shape=[S,N])

    s = inputs

    umprev,upprev = prevstate

    um = (1.-dt/taum)*umprev + (dt/taum)*(gm*tf.matmul(tf.nn.relu(upprev),Jpm*Jpm_mask) + bm)
    up = (1.-dt/taup)*upprev + (dt/taup)*(gp*tf.matmul(tf.nn.relu(upprev),Jpp*Jpp_mask) + bp + tf.matmul(s,wsp))

    return [um,up]

tf.reset_default_graph()

#initial condition
um0_t = tf.zeros(shape=[B,M])
up0_t = tf.placeholder(dtype=np.float32,shape=[B,N],name='up0')

#target
costw_t = tf.placeholder(dtype=tf.float32,shape=[T,B,M],name='costw') #TxBxM
mtarg_t = tf.convert_to_tensor(mtarg.astype(np.float32),name='mtarg') #TxBxM

lr_t = tf.placeholder(dtype=tf.float32,name='lr')
alpha_t = tf.placeholder(dtype=tf.float32,name='alpha')

Jpp0_mask = (Jpp0 > 0).astype(np.float32)
Jpm0_mask = (Jpm0 > 0).astype(np.float32)

Jpp0i = initJpp(Jpp0,types)
Jpm0i = initJpp(Jpm0,types)

#variables
with tf.variable_scope('model'):
    Jpm = tf.get_variable('Jpm',initializer=Jpm0i,trainable=True)
    Jpp = tf.get_variable('Jpp',initializer=Jpp0i,trainable=True)
    bm = tf.get_variable('bm',initializer=0.*np.ones(M,dtype=np.float32))
    bp = tf.get_variable('bp',initializer=0.*np.ones(N,dtype=np.float32))
    taum = tf.get_variable('taum',initializer=tau0*np.ones(M,dtype=np.float32))
    taup = tf.get_variable('taup',initializer=tau0*np.ones(N,dtype=np.float32))
    gm = tf.get_variable('gm',initializer=g0*np.ones(M,dtype=np.float32))
    gp = tf.get_variable('gp',initializer=g0*np.ones(N,dtype=np.float32))
    wsp = tf.get_variable('wsp',initializer=wsp0.astype(np.float32))

    Jpm_mask = tf.get_variable('Jpm_mask',initializer=Jpm0_mask,trainable=False)
    Jpp_mask = tf.get_variable('Jpp_mask',initializer=Jpp0_mask,trainable=False)

#scan
s_t = tf.placeholder(dtype=np.float32,shape=[T,B,S],name='s')

inputs_t = s_t
initial_state_t = [um0_t,up0_t]
scan_out = tf.scan(scan_step,inputs_t,initial_state_t)
um_t,up_t = scan_out
m_t = tf.nn.relu(um_t)
p_t = tf.nn.relu(up_t)

#cost
Jpall = tf.concat([Jpp,Jpm],axis=1)

cost_targ = tf.reduce_sum(tf.pow(costw_t*(m_t-mtarg_t),2))

cost_1827 = alpha1827 * (tf.reduce_sum(tf.abs(p_t[:,0,ind18[0]])) + tf.reduce_sum(tf.abs(p_t[:,0,ind18[1]])) + tf.reduce_sum(tf.abs(p_t[:,1,ind27[0]])) + tf.reduce_sum(tf.abs(p_t[:,1,ind27[1]])))

exec(open("matchnames.py").read()) #generates cost_seg

cost_J = alpha_J * alpha_t * tf.reduce_sum(tf.pow(Jpp*Jpp_mask - Jpp0i,2)) + tf.reduce_sum(tf.pow(Jpm*Jpm_mask - Jpm0i,2))

cost = (cost_targ + cost_1827 + cost_seg + cost_J)/B

train_step = tf.train.RMSPropOptimizer(lr_t).minimize(cost)

Jppmin = -np.infty * np.ones(Jpp.shape)
Jppmax = np.infty * np.ones(Jpp.shape)
Jpmmin = -np.infty * np.ones(Jpm.shape)
Jpmmax = np.infty * np.ones(Jpm.shape)

for qi in range(N):
    if types[qi] == 'exc':
        Jppmin[qi,:] = 0
        Jpmmin[qi,:] = 0
    elif types[qi] == 'inh':
        Jppmax[qi,:] = 0
        Jpmmax[qi,:] = 0

clipJpp_op = tf.assign(Jpp,tf.clip_by_value(Jpp,Jppmin,Jppmax))
clipJpm_op = tf.assign(Jpm,tf.clip_by_value(Jpm,Jpmmin,Jpmmax))
cliptaum_op = tf.assign(taum,tf.clip_by_value(taum,taumin,taumax))
cliptaup_op = tf.assign(taup,tf.clip_by_value(taup,taumin,taumax))
clipgm_op = tf.assign(gm,tf.clip_by_value(gm,0,np.infty))
clipgp_op = tf.assign(gp,tf.clip_by_value(gp,0,np.infty))
clip_op = [clipJpp_op,clipJpm_op,cliptaum_op,cliptaup_op,clipgp_op,clipgm_op]

#initialize session
saver = tf.train.Saver()
try:
    sess
except NameError:
    sess = tf.InteractiveSession()
else:
    sess.close()
    sess = tf.InteractiveSession()
sess.run(tf.global_variables_initializer())
