NEURON {
    SUFFIX CaDynamics
    USEION ca READ ica WRITE cai   
    RANGE tau_Ca, cai_inf, R_Ca_channel, R_removal, cai
}

UNITS {
    (mM) = (milli/liter)
    (nA) = (nanoamp)
    (um3) = (micron3)
    FARADAY = 96500 (coul)
}

PARAMETER {
    f_I = 1e-3 (s/ms)            : Conversion factor for ICa (nA -> nC/ms)
    v = 2000 (um3)               : Effective volume in which [Ca2+]i is changing
    fv = 1e-15 (um3)             : Conversion factor for volume (um3)
    z = 2                         : Valency of calcium ion
    cai_inf = 5e-5 (mM)          : Steady-state calcium concentration
    tau_Ca = 20 (ms)             : Time constant for calcium removal
    fF = 1e6 (nC/mol)            : Conversion factor for Faraday's constant
}

ASSIGNED {
    ica (nA)                     : Calcium current
    R_Ca_channel (mM/ms)         : Rate of Ca2+ entry
    R_removal (mM/ms)            : Rate of Ca2+ removal
}

STATE {
    cai (mM)                     : Intracellular calcium concentration
}

INITIAL {
    cai = cai_inf                : Set steady-state concentration as initial state
}

BREAKPOINT {
    SOLVE state METHOD cnexp
}

DERIVATIVE state {
    : Set rate of Ca2+ entry based on calcium current
    if (ica < 0) {
        R_Ca_channel = - (f_I * ica) / (fv * v * z * fF * FARADAY)
    } else {
        R_Ca_channel = 0
    }

    : Rate of calcium removal
    R_removal = (cai - cai_inf) / tau_Ca

    : Differential equation for calcium concentration
    cai' = R_Ca_channel - R_removal
}
