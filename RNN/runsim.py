#this file is part of zarin_et_al_multilayer_2019
#Copyright (C) 2019 Ashok Litwin-Kumar
#see README for more information

import h5py 
import time
import tensorflow as tf

import matplotlib.pyplot as plt
plt.ion()
plt.rcParams['image.aspect'] = 'auto'

exec(open("setupmodel.py").read())

Nepochspc = max(int(Nepochs/100),1) #show progress every 1% progress interval
track_cost = np.zeros(100)
lastt = time.time()

count = 0
for ei in range(Nepochs):
    #run one step of optimization
    up0 = sigp0*np.abs(np.random.standard_normal([B,N]).astype(np.float32))
    feed_dict = {alpha_t: alpha[ei], s_t: s, costw_t: costw, lr_t: lr[ei], up0_t: up0}
    sess.run(train_step,feed_dict=feed_dict)

    #clip variables
    sess.run(clip_op)

    #show progress
    if (ei % Nepochspc) == 0:
        curt = time.time()
        print("\r" + str(int(100*ei/Nepochs)) + "%,", np.round(curt-lastt,2), "seconds, cost =",track_cost[max(count-1,0)], end="")
        lastt = curt
        track_cost[count],m,p = sess.run([cost_targ,m_t,p_t],feed_dict=feed_dict)

        plt.clf()
        plt.subplot2grid((4,2),(0,0),colspan=2)
        plt.semilogy(track_cost[0:count]/B,color="k")
        plt.ylim(0.1,1.1*np.max(track_cost))
        plt.xlim(0,100)
        plt.ylabel("cost")

        plt.subplot(423)
        plt.imshow(mtarg[:,0,:].T)
        plt.ylabel("MN targets")
        plt.title("FWD")

        plt.subplot(424)
        plt.imshow(mtarg[:,1,:].T)
        plt.title("BWD")

        plt.subplot(425)
        plt.imshow(m[:,0,:].T)
        plt.ylabel("MNs")

        plt.subplot(426)
        plt.imshow(m[:,1,:].T)

        plt.subplot(427)
        plt.imshow(p[:,0,:].T)
        plt.ylabel("PMNs")
        plt.xlabel("timestep")

        plt.subplot(428)
        plt.imshow(p[:,1,:].T)
        plt.xlabel("timestep")

        plt.pause(.0001)
        plt.show()
        plt.tight_layout()

        count += 1
