NEURON {
    SUFFIX Cl_Ca
    USEION ca READ cai       : Reads intracellular calcium concentration [Ca2+]i
    USEION cl WRITE icl CHARGE -1     : Outputs chloride current icl
    RANGE gClCa_max, eClCa, iClCa
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (molar) = (1/liter)			: moles do not appear in units
    (mM)	= (millimolar)
    (nS) = (nanosiemens)
}

PARAMETER {
    v (mV)
    cai (mM)                        : Intracellular calcium concentration [Ca2+]i
    gClCa_max = 4.31e-6 (mho/cm2) 
    :0.0000013 (mho/cm2) : Maximum chloride conductance
    eClCa = -20 (mV)                : Reversal potential for chloride
    KClCa = 1.5e-3 (mM)             : Half-activation concentration
}

ASSIGNED {
    iClCa (mA/cm2)                  : Chloride current
    nClCa                           : Activation variable for chloride channel
    icl (mA/cm2)                    : Chloride current output to NEURON
}

BREAKPOINT {
    nClCa = 1 / (1 + pow((KClCa / cai), 4))
    iClCa = gClCa_max * nClCa * (v - eClCa)
    icl = iClCa                     : Outputs chloride current to NEURON
}
