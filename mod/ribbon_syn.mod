: four-state depressing synapse model for graded membrane potential

NEURON {
    POINT_PROCESS ribbon_syn
    POINTER v_pre
    RANGE e, tau_1A, tau_A3, tau_32, tau_21
    RANGE v_th, v_slp, ca, acm
    RANGE g_max, u, isyn, alpha, beta
    NONSPECIFIC_CURRENT i
}

PARAMETER {
    tau_1A = 2.0 (ms)
:    tau_AI = 10.0 (ms)
:    tau_I3 = 10000 (ms) :added
    tau_A3 = 10000 (ms)
    tau_32 = 2000 (ms)
    tau_21 = 700 (ms)
    P1max = 0.01
    P2max = 0.04
    e = 0 (millvolts)
    g_max = 0.00256 (umho)  : maximal conductance
    v_slp = 10.0 (millvolts) :10
    v_th = -40.0 (millvolts) :-40
    alpha = 1.0 (1/ms)
    beta = 1.1 (1/ms)
}

ASSIGNED {
    v (millivolt)
    v_pre (millivolt)
    i (nanoamp)
    isyn (nanoamp)
    ca
    u
    g
    acm
}

STATE {
    act p1 p2 w
}

BREAKPOINT {
    SOLVE states METHOD euler
:    g = g_max * u*p1
    g = g_max * w
    isyn = g * (v - e)
    i = isyn
}

INITIAL {
 act = 0.047
 p1 = 0.447
 p2 = 0.7
 w = 0.0
:    P3 = 0.9
}

DERIVATIVE states {
    u = (1.0 + tanh((v_pre - v_th)/v_slp))/2.0
    :ca = 0.038/(exp(-(v_pre-(-34))/6) + exp((v_pre-(23))/40))
    :u = ca*ca*ca/(ca*ca*ca + 0.087*0.087*0.087)

    act' = P1max*u*p1/tau_1A - act/tau_A3
:    act' = P1max*u*p1/tau_1A - act/tau_AI
:    inact' = act/tau_AI - inact/tau_I3

:    P3' = act/tau_A3 - P3*(1-p2)*P2max/tau_32  
:    p2' = P3*(1-p2)/tau_32 - p2*(1-p1)*P1max/P2max/tau_21

    p2' = (1-P1max*p1-P2max*p2-act)*(1-p2)/tau_32 - p2*(1-p1)*P1max/P2max/tau_21
    p1' = p2*(1-p1)/tau_21 - u*p1/tau_1A

    w' = alpha*u*p1 -beta*w
    acm = acm + alpha*u*p1
}