INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX ikxrod
    USEION k READ ek WRITE ik
    :NONSPECIFIC_CURRENT ikx
    RANGE gkxbar
    RANGE ninf_kx, taun_kx
    RANGE ikx
}

PARAMETER {
    gkxbar = 2.82e-6 (mho/cm2) 
    :0.00000104 (mho/cm2) : Maximum conductance
    ekx = -80 (mV)           : : Reversal potential
    v (mV)
    tau_kx (ms)        : Time constant
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (milli/liter)
    (nS) = (nanosiemens)
} 

ASSIGNED {
    ik (mA/cm2)
    ikx (mA/cm2)
    ek (mV)                    : Reversal potential
    ninf_kx   
    taun_kx (ms)
    a (1/ms)
    b (1/ms)
    :a0 (1/ms)
    :b0 (1/ms)
}

STATE {
    n_kx
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ikx = gkxbar * n_kx * (v - ek)
    ik = ikx
} 

INITIAL {
    trates_kx(v)
    n_kx = ninf_kx
}

DERIVATIVE states {   
    trates_kx(v)      
    n_kx' = a * (1 - n_kx) - b * n_kx
}

PROCEDURE trates_kx(vm) { 
    LOCAL a, b
    a = alpha_kx(vm, -50, 0.00066, 5.7)
    b = beta_kx(vm, -50, 0.00066, 5.7)
    ninf_kx = a / (a + b)
    taun_kx = 1 / (a + b)

}

FUNCTION alpha_kx(v, vhalf, a0, s) { 
    alpha_kx = a0 * (exp((v - vhalf) / 2*s))
}   

FUNCTION beta_kx(v, vhalf, b0, s) { 
    beta_kx = b0 * (exp(-(v - vhalf) / 2*s))
}