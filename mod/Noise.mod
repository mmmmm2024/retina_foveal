NEURON {
    POINT_PROCESS Noise
    RANGE i,del,dur,f0,f1,r,torn,std,bias
    ELECTRODE_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    del=0    (ms)
    dur=5000   (ms)
    torn=500  (ms)
    std=0.2   (nA)
    f0=0.2    (nA)
    f1=0.8    (nA)
    r =60
    bias = 0 (nA)
}

ASSIGNED {
    ival (nA)
    i (nA)
    noise (nA)
:amp (nA)
}

INITIAL {
    i = 0
}

PROCEDURE seed(x) {
    set_seed(x)
}

BEFORE BREAKPOINT {
    noise = normrand(0,std*1(/nA))*1(nA)
:amp = f0 + 0.5*(f1-f0)*(tanh((t-torn)/(r/3)/(1(ms))-3)+1)
:ival = amp + noise + bias
    ival = noise + bias
}

BREAKPOINT {
    i = ival
}