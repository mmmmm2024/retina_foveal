:Yonemoto /* Kca channel*/
NEURON {
    SUFFIX IKCayone
    USEION ca READ cai           
    USEION k WRITE ik            
    RANGE gKca, EKca, ik         
    RANGE mKca                   
    RANGE Cahalf_mKca            
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (molar) = (1/liter)			: moles do not appear in units
    (mM)	= (millimolar)
}

PARAMETER {
    gKca = 0.005 (mho/cm2)           
    EKca = -80 (mV)               
    Cahalf_mKca = 0.00032 (mM)    
}

ASSIGNED {
    v (mV)                         
    cai (mM)                       
    ik (mA/cm2)                    
    mKca                           
}

BREAKPOINT {
    mKca = mKca_RPR(cai)           
    ik = gKca * mKca * mKca *mKca * mKca *(v - EKca) 
}

INITIAL {
    mKca = 0.5
}

FUNCTION mKca_RPR(ca) {
    mKca_RPR = 1.0 / (1.0 + (Cahalf_mKca / ca) * (Cahalf_mKca / ca) * (Cahalf_mKca / ca) * (Cahalf_mKca / ca))
}
