TITLE HH k channel channel
: Hodgkin - Huxley A channel


NEURON {
	SUFFIX A
	USEION k READ ek WRITE ik
	RANGE gIAbar, ik
	GLOBAL minf, hinf, cfunc, h1tau, h2tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gIAbar=.00026 (mho/cm2) <0,1e9>
	ek = -77 (mV) :suggested, default set by NEURON
}

STATE {
	m h1 h2 c
}

ASSIGNED {
	v (mV)
	celsius (degC) : 16
	ik (mA/cm2)
	minf
    hinf
    h1tau (ms)
    h2tau (ms)
    cfunc (mV)
}

INITIAL {
    rate(v*1(/mV))
	m = minf
    h1 = hinf
    h2 = hinf
    c = cfunc
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    ik = gIAbar*m*(c*h1+(1-c)*h2)*(v - ek)
}

DERIVATIVE states {	: exact when v held constant
	rate(v*1(/mV))
	m' = (minf - m)/1
    h1' = (hinf - h1)/h1tau
    h2' = (hinf - h2)/h2tau
}

UNITSOFF
FUNCTION m_inf(v(mV)) {
    m_inf = 1/(1+exp(-( v + 10 ) / 7))
}

FUNCTION h_inf(v(mV)) {
    h_inf = 0.83 * (1/(1+exp(( v + 40.5 ) / 2))) + 0.17
}

FUNCTION h1_tau(v(mV)) { :phi=0.039
    h1_tau = 1/(25 - 20/(1 + exp(-( v + 35 ) / 6)))
}

LOCAL tmp1, tmp2

FUNCTION h2_tau(v(mV)) { :phi=0.039
    tmp1 = (((v+17)*(v+17)) / 4) + 26 
    tmp2 = 100
    if(tmp1<tmp2){
        h2_tau = tmp1
    }else{
        h2_tau = tmp2
    }
}

FUNCTION c_func(v(mV)) {
    c = 1/(1+exp(-( v + 45 ) / 15))
}

PROCEDURE rate(v(mV)) {
	TABLE minf,hinf,h1tau,h2tau,c FROM -100 TO 100 WITH 200
		h1tau = h1_tau(v)
        h2tau = h2_tau(v)
        minf = m_inf(v)
        hinf = h_inf(v)
        cfunc = c_func(v)
}
UNITSON
