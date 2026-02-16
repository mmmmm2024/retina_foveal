TITLE L-type calcium channel for Tiger Salamander Bipolar cell
:
: Modified from Fohlmeister et al, 1990, Brain Res 510, 343-345
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
 SUFFIX cabip
 USEION ca READ cai, eca, cao WRITE ica
 RANGE gcabar
 RANGE m_inf, tau_m, m_exp
}


UNITS {
 (molar) = (1/liter)
 (mM) = (millimolar)
 (mA) = (milliamp)
 (mV) = (millivolt)
}

PARAMETER {
 gcabar	= 0.00005 (mho/cm2) : 0.00055
 eca (mV)
 cao = 2.5 (mM) :1.8
 cai = 0.0001 (mM)
 dt (ms)
 v  (mV)

}

STATE {
 m h
}

INITIAL {
 m = 0.24
 h = 1.0
}

ASSIGNED {
 ica (mA/cm2)
 m_inf
 tau_m
 m_exp
}

BREAKPOINT {
 SOLVE states
 ica = gcabar * m*m*m*m *h* (v - eca)
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
	
 am = 12000*(120-v) / ((exp(-(v-120)/25)) - 1)
 bm = 40000 / (exp((v+68)/25) + 1)
 h = exp(-(v-50)/11)/(exp(-(v-50)/11) + 1)

 tau_m = 1 / (am + bm)
 m_inf = am * tau_m

: State vars to inifinity
 m_exp = 1 - exp(-dt/tau_m)

}

UNITSON
