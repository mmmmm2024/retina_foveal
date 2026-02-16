:Yonemoto /* Kv channel */
NEURON {
    SUFFIX IKvyone
    USEION k WRITE ik     : Uses potassium ion and writes potassium current
    RANGE gKv, EKv, ik            : Parameters for maximum conductance, reversal potential, and current
    RANGE mKv, m_inf, tau_m       : Activation variable, steady-state value, and time constant for m
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mS) = (millisiemens)
}

PARAMETER {
    gKv = 0.01 (mho/cm2)           : Maximum conductance of the Kv channel
    EKv = -80 (mV)                : Reversal potential for potassium
}

ASSIGNED {
    v (mV)                        : Membrane potential
    ek (mV)                       : Reversal potential for potassium
    ik (mA/cm2)                   : Potassium current
    m_inf                         : Steady-state activation variable
    tau_m (ms)                    : Time constant for m
}

STATE {
    mKv                           : Activation variable for the Kv channel
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gKv * mKv * mKv *  mKv * mKv * (v - EKv) : Voltage-dependent potassium current
}

DERIVATIVE states {
    rates(v)
    mKv' = (m_inf - mKv) / tau_m  : Differential equation for mKv
}

INITIAL {
    rates(v)
    mKv = m_inf                   : Set initial value of mKv to its steady state
}

PROCEDURE rates(v (mV)) {
    LOCAL alpha_m, beta_m
    alpha_m = alpha_mKv_RPR(v)
    beta_m = beta_mKv_RPR(v)
    m_inf = alpha_m / (alpha_m + beta_m)
    tau_m = 1.0 / (alpha_m + beta_m)
}

FUNCTION alpha_mKv_RPR(v (mV)) (/ms) {
    alpha_mKv_RPR = 0.001 * 5 * (20.0 - v) / (exp((20.0 - v) / 22.0) - 1.0)
}

FUNCTION beta_mKv_RPR(v (mV)) (/ms) {
    beta_mKv_RPR = exp(-v / 80.0) / 16.0
}
