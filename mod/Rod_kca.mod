NEURON {
    SUFFIX KCa
    USEION ca READ cai       
    USEION k WRITE ik        : Kalium current
    RANGE gKCa_max, eKCa, iKCa
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (milli/liter)
    (nS) = (nanosiemens)
}

PARAMETER {
    v (mV)
    cai (mM)                        : [Ca2+]i
    gKCa_max = 1.65e-5 (mho/cm2)  
    :0.000005 (mho/cm2)   : Max Kalium conductance
    eKCa = -80 (mV)                 
    KKCa = 3.2e-4 (mM)              :half-activation concentration
}

ASSIGNED {
    iKCa (mA/cm2)                   
    nKCa                            
    ik (mA/cm2)                     
}

BREAKPOINT {
    nKCa = 1 / (1 + pow((KKCa / cai), 4))
    iKCa = gKCa_max * nKCa^4 * (v - eKCa)
    ik = iKCa                       
}
