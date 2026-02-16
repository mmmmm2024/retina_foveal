TITLE L-type calcium channel for Tiger Salamander Bipolar cell
:
: Modified from Fohlmeister et al, 1990, Brain Res 510, 343-345
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
 SUFFIX kcabip
 USEION k READ ek WRITE ik
 USEION ca READ cai
 RANGE gkcabar
 RANGE m_inf, tau_m, m_exp
 RANGE ikca
}


UNITS {
 (molar) = (1/liter)
 (mM) = (millimolar)
 (mA) = (milliamp)
 (mV) = (millivolt)
}

PARAMETER {
 gkcabar = 0.00425
 cai (mM)
 eca (mV)
 ek  (mV)
 dt  (ms)
 v   (mV)
}

STATE {
 m
}

INITIAL {
 m = 0.45
}

ASSIGNED {
 ik   (mA/cm2)
 ikca (mA/cm2)
 m_inf
 tau_m
 m_exp
}

BREAKPOINT {
 SOLVE states
 ikca = gkcabar * m*m*((cai)/(cai+0.2)) * (v - ek)
 ik = ikca
}

PROCEDURE states() {	: exact when v held constant
 evaluate_fct(v)
 m = m + m_exp * (m_inf - m)

 VERBATIM
 return 0;
 ENDVERBATIM
}

UNITSOFF

PROCEDURE evaluate_fct(v(mV)) { LOCAL am, bm

 am = (100* (230-v)) / ((exp((230-v)/52)) - 1)
 bm = 120 * (exp((-v/95)))

 tau_m = 1 / (am + bm)
 m_inf = am * tau_m

: State vars to inifinity
 m_exp = 1 - exp(-dt/tau_m)
}

UNITSON
