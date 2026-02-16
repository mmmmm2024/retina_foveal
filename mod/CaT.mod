TITLE low-voltage-activated (T-type) calcium channel for STN neuron

COMMENT
 modeled by Otsuka et al., 2004
 implemented in NEURON by Kitano, 2017
ENDCOMMENT

UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (mM) = (milli/liter)
}

NEURON {
 SUFFIX CaT
 USEION ca READ cai,cao WRITE ica
 RANGE gmax, iCaT
}

PARAMETER {
 v (mV)
 dt (ms)
 cai (mM)
 cao (mM)
 gmax  = 0.001 (mho/cm2)
 iCaT  = 0.0 (mA/cm2)
 e = 130 (mV)
}

STATE {
 m h
}

ASSIGNED { 
 ica (mA/cm2)
 minf
 taum
 hinf
 tauh
}

BREAKPOINT {
 SOLVE states METHOD cnexp
 ica  = gmax*m*m*h*(v-e)
 iCaT = ica
}

UNITSOFF

INITIAL {
 settables(v)
 m = minf
 h = hinf
}

DERIVATIVE states {  
 settables(v)
 m' = (minf - m)/taum
 h' = (hinf - h)/tauh
}

PROCEDURE settables(v) {
        TABLE minf, taum, hinf, tauh FROM -100 TO 100 WITH 400

	minf = 1.0 / (1.0 + exp(-(v + 56.0)/6.7))
	taum = 5.0 + 0.33/(exp((v + 27.0)/10.0) + exp(-(v + 102.0)/15.0))
	hinf = 1.0 / (1.0 + exp((v + 85.0)/5.8))
	tauh = 400.0/(exp((v + 50.0)/15.0) + exp(-(v + 50.0)/16.0))
}

UNITSON
