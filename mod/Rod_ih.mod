INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX ihrod
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT ih
    RANGE ghbar
    RANGE ninf_h, taun_h
    RANGE ih
}

PARAMETER {
    ghbar = 8.27e-6 (mho/cm2) 
    :0.0000025 (mho/cm2) : Maximum conductance
    eh = -32 (mV)           : : Reversal potential
    v (mV)
    :tau_kv = 100 (ms)        : Time constant
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (milli/liter)
    (nS) = (nanosiemens)
}

ASSIGNED {
    ik (mA/cm2)
    ih (mA/cm2)
    ek (mV)                    : Reversal potential
    ninf_h   
    taun_h (ms)
    a (1/ms)
    b (1/ms)
    :a0 (1/ms)
    :b0 (1/ms)
}

STATE {
    n_h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ih = ghbar * n_h *(v - ek)
    ik = ih
} 

INITIAL {
    trates_h(v)
    n_h = ninf_h

}

DERIVATIVE states {   
    trates_h(v)      
    n_h' = a * (1 - n_h) - b * n_h
}

PROCEDURE trates_h(vm) { 
    LOCAL a, b
    a = alpha_h(vm, -75, -5.3)
    b = beta_h(vm, -75, -5.3)
    ninf_h = a / (a + b)
    taun_h = 1 / (a + b)
}

FUNCTION alpha_h(v, vhalf, s) { 
    alpha_h = 0.001 * (exp((v - vhalf) / 2*s))
}   

FUNCTION beta_h(v, vhalf, s) { 
    beta_h = 0.001 * (exp(-(v - vhalf) / 2*s))
}