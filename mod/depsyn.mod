: three-state depressing synapse model for graded membrane potential

NEURON {
    POINT_PROCESS depsyn
    POINTER v_pre
    RANGE e, tau_e, tau_r, v_th, v_slp, s, g_max, u
    NONSPECIFIC_CURRENT i
}

PARAMETER {
    tau_e = 2.7 (ms)
    tau_r = 800 (ms)
    e = -70 (millvolts)
    g_max = 0.00256 (umho)  : maximal conductance
    v_slp = 20.0 (millvolts) :10
    v_th = -40.0 (millvolts) :-40
    u = 0.4 :0.005
}

ASSIGNED {
    v (millivolt)
    v_pre (millivolt)
    i (nanoamp)
    s
    g
}

STATE {
    eff rec
}

BREAKPOINT {
    SOLVE states METHOD euler
    g = g_max * eff
    i = g * (v - e)
}

INITIAL {
    eff = 0.01
    rec = 1.0
}

DERIVATIVE states {
    s = (1.0+tanh((v_pre - v_th)/v_slp))/2.0
    eff' = -eff/tau_e + u*s*rec
    rec' = (1.0-rec-eff)/tau_r - u*s*rec
}
