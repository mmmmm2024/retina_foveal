TITLE HH style channels for spiking retinal ganglion cells
:
: Modified from Fohlmeister et al, 1990, Brain Res 510, 343-345
: by TJ Velte March 17, 1995
: must be used with calcium pump mechanism, i.e. capump.mod
:
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX spike2
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	USEION ca READ cai, eca, cao WRITE ica
	RANGE gnabar, gkbar, gcabar, gkcbar, gkc
	RANGE m_inf, h_inf, n_inf, c_inf
	RANGE tau_m, tau_h, tau_n, tau_c
	RANGE m_exp, h_exp, n_exp, c_exp
	RANGE idrk, icak
	GLOBAL ca50
}


UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar	= 0.08	(mho/cm2)
	gkbar	= 0.05 (mho/cm2)
	gcabar	= 0.0012	(mho/cm2)
	gkcbar	= 0.00005 (mho/cm2)
	ena	= 65	(mV)
	ek	= -100	(mV)
	eca		(mV)
	cao	= 1.8	(mM)
	cai     = 0.0001 (mM)
	dt              (ms)
	v               (mV)
	ca50    = 0.001	(mM)
}

STATE {
	m h n c
}

INITIAL {
: The initial values were determined at a resting value of -66.3232 mV in a single-compartment
:	m = 0.0155
:	h = 0.9399
:	n = 0.0768
:	p = 0.0398
:	q = 0.4526
:	c = 0.0016
: at -60 mV
        m = 0.0345
        h = 0.8594
        n = 0.1213
        c = 0.0038
}

ASSIGNED {
	ina	(mA/cm2)
	ik	(mA/cm2)
         idrk    (mA/cm2)
         icak    (mA/cm2)
		 ikm	 (mA/cm2)
	ica	(mA/cm2)
	m_inf h_inf n_inf c_inf
	tau_m tau_h tau_n tau_c
	m_exp h_exp n_exp c_exp
	gkc

}

BREAKPOINT {
	SOLVE states
	ina = gnabar * m*m*m*h * (v - ena)
	idrk = gkbar * n*n*n*n * (v - ek)
:	gkc=((cai / ca50)/ (1 + (cai / ca50)))
	gkc=(cai*cai)/ (ca50*ca50 + cai*cai)
	icak = gkcbar * gkc * (v - ek)
	ik = idrk + icak
	ica = gcabar * c*c*c * (v - eca)

}

PROCEDURE states() {	: exact when v held constant
	LOCAL sigmas
	evaluate_fct_nak(v)
	evaluate_fct_ca(v)
	m = m + m_exp * (m_inf - m)
	h = h + h_exp * (h_inf - h)
	n = n + n_exp * (n_inf - n)
	c = c + c_exp * (c_inf - c)

	VERBATIM
	return 0;
	ENDVERBATIM

}

UNITSOFF
PROCEDURE evaluate_fct_ca(v(mV)) { LOCAL a,b
:CA channel
	a = (-1.362 * (v+13)) / (exp(-0.1*(v+13)) - 1)
	b = 45.41 * exp(-1*(v + 38)/18)
	tau_c = 1 / (a + b)
	c_inf = a * tau_c
	c_exp = 1 - exp(-dt/tau_c)	

}
PROCEDURE evaluate_fct_nak(v(mV)) { LOCAL a,b
	
:NA m
	:a = (-2.725 * (v+13)) / (exp(-0.1*(v+35)) - 1)
	a = (-2.725 * (v+35)) / (exp(-0.1*(v+35)) - 1)
	b = 90.83 * exp(-1*(v+60)/18)
	tau_m = 1 / (a + b)
	m_inf = a * tau_m
	m_exp = 1 - exp(-dt/tau_m)

:NA h
	a = 1.817 * exp(-1*(v+52)/20)
	b = 27.25 / ( 1 + exp(-0.1 *(v+22)))
	tau_h = 1 / (a + b)
	h_inf = a * tau_h
	h_exp = 1 - exp(-dt/tau_h)

:K n (non-inactivating, delayed rectifier)
	a = (-0.09575 * (v+37)) / (exp(-0.1*(v+37)) - 1)
	b = 1.915 * exp(-1*(v + 47)/80)
	tau_n = 1 / (a + b)
	n_inf = a * tau_n
	n_exp = 1 - exp(-dt/tau_n)

}

UNITSON
