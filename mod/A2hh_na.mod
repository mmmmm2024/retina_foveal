TITLE HH sodium channel
: Hodgkin - Huxley squid sodium channel
: for Riecke 2014

NEURON {
	SUFFIX HHna
	USEION na READ ena WRITE ina
	RANGE gnabar, ina
	GLOBAL minf, hinf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar=0.0026 (mho/cm2) <0,1e9>
    ena = 50(mv)
}

STATE {
	m h
}

ASSIGNED {
	v (mV)
	celsius (degC) : 6.3
	ina (mA/cm2)
	minf hinf
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    ina = gnabar*m*m*m*h*(v - ena)
}

DERIVATIVE states {
	rates(v)
    m' = (minf - m)/0.01
	h' = (hinf - h)/0.5
}

FUNCTION m_inf(v(mV)) {
    m_inf = 1/(1+exp(-(v + 48) / 5))
}

FUNCTION h_inf(v(mV)) {
    h_inf = 1/(1+exp((v + 49.5) / 2))
}


PROCEDURE rates(v(mV)) {
	TABLE minf, hinf FROM -100 TO 100 WITH 200
	minf = m_inf(v)
    hinf = h_inf(v)
}
