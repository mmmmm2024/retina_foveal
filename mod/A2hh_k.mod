TITLE HH k channel channel
: Hodgkin - Huxley M-type k channel


NEURON {
	SUFFIX Ik
	USEION k READ ek WRITE ik
	RANGE gIkbar, ik
	GLOBAL minf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gIkbar=.00026 (mho/cm2) <0,1e9>
	ek = -77(mV) :suggested, default set by NEURON
}

STATE {
	m
}

ASSIGNED {
	v (mV)
	celsius (degC) : 16
	ik (mA/cm2)
	minf
}

INITIAL {
    rate(v*1(/mV))
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    ik = gIkbar*m*(v - ek)
}

DERIVATIVE states {	: exact when v held constant
	rate(v*1(/mV))
	m' = (minf - m)/50
}

UNITSOFF
FUNCTION m_inf(v(mV)) {
    m_inf = 1/(1+exp(-(v + 40 ) / 4))
}

PROCEDURE rate(v(mV)) {
	TABLE minf FROM -100 TO 100 WITH 200
        minf = m_inf(v)
}
UNITSON
