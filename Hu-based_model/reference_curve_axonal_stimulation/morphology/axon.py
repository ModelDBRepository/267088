def makeAxon(h, N_ranvier=3, Return_Connected=False):
    # from neuron import h
    import numpy as np
    from myelin import makeMyelin
    from ranvier import makeRanvier
    N_internodes = N_ranvier+1

    axon = []
    nodes_of_ranvier = []
    internodes = []
    for i in np.arange(N_internodes):
        m = makeMyelin(h, i)
        axon.append(m)
        internodes.append(m)
        if i < N_ranvier:
            r = makeRanvier(h, i)
            axon.append(r)
            nodes_of_ranvier.append(r)

    if Return_Connected == True:
        for i in np.arange(len(axon)-1):
            axon[i+1].connect(axon[i](1), 0)

    return axon, nodes_of_ranvier, internodes 