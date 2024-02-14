


def cosine_taper(left_section, tapering_section, right_section):
    import numpy as np

    # tapering_section.nseg = 31 #<---Must be odd number
    if tapering_section.nseg < int(tapering_section.L):
        tapering_section.nseg = 11 + 2*int(tapering_section.L) #<---Must be odd number
    
    diam_left = left_section(1).diam
    diam_right = right_section(0).diam
    
    if diam_left != diam_right:
        if diam_left > diam_right:
            for seg in tapering_section.allseg():
                seg.diam = diam_right + (diam_left - diam_right)*(0.5*(1.0+np.cos(np.pi*seg.x)))

        elif diam_right > diam_left:
            for seg in tapering_section.allseg():
                seg.diam = diam_left + (diam_right - diam_left)*(0.5*(1.0+np.cos(np.pi*(1.0 + seg.x))))

    elif diam_left == diam_right:
        tapering_section.diam = diam_left


def conical_taper(left_section, tapering_section, right_section):

    # tapering_section.nseg = 31 #<---Must be odd number
    if tapering_section.nseg < int(tapering_section.L):
        tapering_section.nseg = 11 + 2*int(tapering_section.L) #<---Must be odd number
    
    diam_left = left_section(1).diam
    diam_right = right_section(0).diam
    
    if diam_left != diam_right:
        for seg in tapering_section.allseg():
            # seg.diam = diam_right + (diam_left - diam_right)*(1.0 - seg.x)
            seg.diam = diam_left + (diam_right - diam_left)*(seg.x)


    elif diam_left == diam_right:
        tapering_section.diam = diam_left


def conical_taper_create_Section(left_Section, right_Section, newSection_length=11.0, name=None):
    from neuron import h
    import numpy as np

    if name == None:
        print('\n\n *********ERROR: NO TAPERING PERFORMED. Must provide name for new tapered Section\n\n')
        exit()
    else:
        sec = h.Section(name=name)
        
        sec.L = newSection_length
        sec.nseg = 21 #<---Must be odd number
        if sec.nseg < int(newSection_length):
            sec.nseg = 1 + 2*int(newSection_length) #<---Must be odd number
        
        diam_left = left_Section(1).diam
        diam_right = right_Section(0).diam
        
        if diam_left != diam_right:
            if diam_left > diam_right:
                for seg in sec.allseg():
                    seg.diam = diam_right + (diam_left - diam_right)*(1.0 - seg.x)

            elif diam_right > diam_left:
                for seg in sec.allseg():
                    seg.diam = diam_right + (diam_left - diam_right)*(1.0 - seg.x)


        elif diam_left == diam_right:
            sec.diam = diam_left

        sec.connect(left_section(1), 0)
        right_Section.connect(sec(1), 0)
        return newSection


def taper𝛽_Connect𝛼𝛽(𝛽, 𝛼):
    import numpy as np
    if 𝛽.name()=='dend' and 𝛼.name()=='soma':
        for k in 𝛽.allseg():
            lambda_𝛽 = 0.01
            # k.diam = 𝛼.diam + (𝛽.diam - 𝛼.diam)*0.5*(1 + np.cos(np.pi*k.x)**2)
            k.diam = 𝛽.diam + (𝛼.diam - 𝛽.diam)/(1.0 + np.exp((0.95-k.x)/lambda_𝛽))
        𝛼.connect(𝛽(1), 0)
        return 𝛽,𝛼 

    elif 𝛽.name()=='hill' and 𝛼.name()=='soma':
        for k in 𝛽.allseg():
            lambda_𝛽 = 0.01
            # k.diam = 𝛼.diam + (𝛽.diam - 𝛼.diam)*0.5*(1 + np.cos(np.pi*k.x)**2)
            k.diam = 𝛽.diam + (𝛼.diam - 𝛽.diam)/(1.0 + np.exp((0.95-k.x)/lambda_𝛽))
        𝛽.connect(𝛼(1), 0)
        return 𝛽,𝛼 

    else:
        print('\n\n *********ERROR: NO TAPERING PERFORMED!!! \n\n')


        