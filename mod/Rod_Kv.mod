INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX ikvrod
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT ikv
    RANGE gkvbar
    RANGE ninf_kv, taun_kv
    RANGE ikv
}

PARAMETER {
    gkvbar = 3.31e-5 (mho/cm2)  
    :0.00001 (mho/cm2) : Maximum conductance
    ekv = -80 (mV)           : : Reversal potential
    v (mV)
    tau_kv = 100 (ms)        : Time constant
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (milli/liter)
    (nS) = (nanosiemens)
}

ASSIGNED {
    ik (mA/cm2)
    ikv (mA/cm2)
    ek (mV)                    : Reversal potential
    ninf_kv   
    taun_kv (ms)
    a (1/ms)
    b (1/ms)
    :a0 (1/ms)
    :b0 (1/ms)
}

STATE {
    n_kv
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ikv = gkvbar * n_kv * n_kv * n_kv * n_kv * (v - ek)
    ik = ikv
} 

INITIAL {
    trates_kv(v)
    n_kv = ninf_kv
}

DERIVATIVE states {   
    trates_kv(v)      
    n_kv' = a * (1 - n_kv) - b * n_kv
}

PROCEDURE trates_kv(vm) { 
    LOCAL a, b
    a = alpha_kv(vm, 0.005)
    b = beta_kv(vm, 0.0625)
    ninf_kv = a / (a + b)
    taun_kv = 1 / (a + b)

}

FUNCTION alpha_kv(v, a0) { 
    alpha_kv = a0 * (exp(22)-1) :(20-v)/((20-v)/(exp(22)-1))
}

FUNCTION beta_kv(v, b0) { 
    beta_kv = b0 * (exp(- v / 80))
}