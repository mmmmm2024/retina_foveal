NEURON {
    POINT_PROCESS ribbon_syn_R2RB
    POINTER v_pre
    RANGE e, tau_1A, tau_A3, tau_32, tau_21
    RANGE v_th, v_slp, ca, acm
    RANGE g_max, u, isyn, alpha, beta
    RANGE P1max, P2max, i_on_threshold
    NONSPECIFIC_CURRENT i
}

PARAMETER {
    tau_1A = 2.0 (ms)
    tau_A3 = 10000 (ms)
    tau_32 = 2000  (ms)
    tau_21 = 700   (ms)

    P1max  = 0.01
    P2max  = 0.04

    g_max = 0.00256 (umho)
    e     = 0 (millvolts)

    v_slp = 10.0 (millvolts)
    v_th  = -40.0 (millvolts)

    alpha = 1.0 (1/ms)
    beta  = 1.1 (1/ms)
 
    :i_on_threshold = -0.01 (nA)
}

ASSIGNED {
    v (millvolts)           
    v_pre (millvolts)       
    i (nanoamp)          
    isyn (nanoamp)         
    ca                
    u                 
    g           
    acm              
    :u_light           
}

STATE {
    act          
    p1               
    p2                
    w                
}

INITIAL {
    act = 0.0
    p1  = 0.015
    p2  = 0.58
    w   = 0.0
    acm = 0
}

BREAKPOINT {
    SOLVE states METHOD euler

    : compute conductance and current
    g    = g_max * w
    isyn = g * (v - e)
    i    = isyn
}

DERIVATIVE states {
    u = (1.0 + tanh((v_pre - v_th)/v_slp)) / 2.0
    act' = P1max*u*p1/tau_1A - act/tau_A3
    p2' = (1-P1max*p1-P2max*p2-act)*(1-p2)/tau_32 - p2*(1-p1)*P1max/P2max/tau_21
    p1' = p2*(1-p1)/tau_21 - u*p1/tau_1A

    w' = alpha*(1-u)*p1 - beta*w
    acm = acm + alpha*(1-u)*p1

    :act' = P1max*u*p1/tau_1A - act/tau_A3

    :p2' = (1 - P1max*p1 - P2max*p2 - act) * (1 - p2)/tau_32 - p2*(1 - p1) * P1max/P2max / tau_21

    :p1' = p2*(1 - p1)/tau_21 - u*p1/tau_1A

    :w' =  alpha * (1 - w) - beta * (1-u) * p1 * w

    :acm = acm + alpha * (1-u) * p1 * w
}