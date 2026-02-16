:yonemoto /* Kx channel */
NEURON {
    SUFFIX IKxyone
    USEION k READ ek WRITE ik        : Uses potassium ion and writes potassium current
    RANGE gKx, EKx, ik               : Parameters for maximum conductance, reversal potential, and current
    RANGE nKx, n_inf, tau_n          : Activation variable, steady-state value, and time constant for nKx
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mS) = (millisiemens)
}

PARAMETER {
    gKx = 0.00085 (mho/cm2)              : Maximum conductance of the Kx channel
    EKx = -80 (mV)                   : Reversal potential for Kx current
    ao_nKx = 0.00066 (/ms)               : Rate constant for Kx activation
    Vhalf_nKx = -49.9 (mV)           : Half-activation voltage
    S_nKx = 5.7 (mV)                 : Slope factor for activation
}

ASSIGNED {
    v (mV)                           : Membrane potential
    ek (mV)                          : Reversal potential for potassium
    ik (mA/cm2)                      : Potassium current
    n_inf                            : Steady-state activation variable
    tau_n (ms)                       : Time constant for nKx
}

STATE {
    nKx                               : Activation variable for the Kx channel
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gKx * nKx * (v - EKx)      : Calculate Kx current based on activation variable
}

DERIVATIVE states {
    rates(v)
    nKx' = (n_inf - nKx) / tau_n     : Differential equation for nKx
}

INITIAL {
    rates(v)
    nKx = n_inf                      : Set initial value of nKx to its steady state
}

PROCEDURE rates(v (mV)) {
    LOCAL alpha_n, beta_n
    alpha_n = alpha_nKx_RPR(v)
    beta_n = beta_nKx_RPR(v)
    n_inf = alpha_n / (alpha_n + beta_n)
    tau_n = 1.0 / (alpha_n + beta_n)
}

FUNCTION alpha_nKx_RPR(v (mV)) (/ms) {
    alpha_nKx_RPR = 0.001 * ao_nKx * exp((v - Vhalf_nKx) / (2 * S_nKx))
}

FUNCTION beta_nKx_RPR(v (mV)) (/ms) {
    beta_nKx_RPR = 0.001 * ao_nKx * exp(-(v - Vhalf_nKx) / (2 * S_nKx))
}
