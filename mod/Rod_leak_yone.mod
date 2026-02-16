NEURON {
    SUFFIX ILeakyone
    NONSPECIFIC_CURRENT il         : Defines a non-specific leak current
    RANGE glbar, El             : Parameters for conductance, reversal potential, and current
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mS) = (millimho) :(millisiemens)
}

PARAMETER {
    glbar = 0.52 (mS/cm2)             : Leak conductance
    El = -74 (mV)                  : Reversal potential for leak current
}

ASSIGNED {
    v (mV)
    gl (mho/cm2)                         : Membrane potential
    il (mA/cm2)                     : Leak current
}

BREAKPOINT {
    gl = 1e-3 * glbar
    il = gl * (v - El)              : Calculate the leak current based on membrane potential
}
