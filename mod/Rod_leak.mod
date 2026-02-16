NEURON {
    SUFFIX leak
    NONSPECIFIC_CURRENT i         : leak current
    RANGE g_leak, e_leak, i       
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (milli/liter)
    (nS) = (nanosiemens)
}

PARAMETER {
    v (mV)                         : membrane potential
    g_leak = 0.52e-9 (mho) 
    :0.00000052 (mho/cm2)  : conductance
    e_leak = -74 (mV)              : reverse potential
}

ASSIGNED {
    i (nA)                         : leak current
}

BREAKPOINT {
    i = g_leak * (v - e_leak)      
}
