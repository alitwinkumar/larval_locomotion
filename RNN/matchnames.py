#this file is part of zarin_et_al_multilayer_2019
#Copyright (C) 2019 Ashok Litwin-Kumar
#see README for more information

matchpair1 = []
matchpair2 = []
for ii in range(N):
    pname = pnames[ii]
    endname = pname[-2:]
    startname = pname[:-2]
    matchname = ""
    if endname == "t3":
        matchname = startname+"a1"
    elif endname == "a1":
        matchname = startname+"a2"
    elif endname == "a2":
        matchname = startname+"a3"
    elif endname == "a3":
        matchname = startname+"a4"
    else:
        print("error, invalid name:",pname)

    if matchname in pnames:
        jj = np.where(matchname == pnames)[0][0]
        matchpair1.append(ii)
        matchpair2.append(jj)

Npairs = len(matchpair1)

matchpairs = np.vstack([np.array(matchpair1),np.array(matchpair2)]).T

t10 = np.min(seg1inds)
t1f = np.max(seg1inds)
t20 = np.min(seg2inds)
t2f = np.max(seg2inds)
cost_seg = 0.
for qq in range(Npairs):
    i1 = matchpairs[qq,0] #index of A1 neuron
    i2 = matchpairs[qq,1] #index of A2 neuron
    cost_seg += tf.reduce_sum(tf.pow(p_t[t10:t1f,0,i2] - p_t[t20:t2f,0,i1],2)) #A2 goes first, then A1
    cost_seg += tf.reduce_sum(tf.pow(p_t[t10:t1f,1,i1] - p_t[t20:t2f,1,i2],2))

cost_seg = alpha_seg*alpha_t*cost_seg


