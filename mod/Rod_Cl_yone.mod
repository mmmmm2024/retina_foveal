:Yonemoto /* Cl channel */
NEURON {
    SUFFIX IClyone
    USEION ca READ cai          
    USEION cl WRITE icl :CHARGE -1
    RANGE gCl, ECl, icl        
    RANGE mCl                   
    RANGE Cahalf_mCl            : Half-activation concentration for mCl
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (molar) = (1/liter)			: moles do not appear in units
    (mM)	= (millimolar)
}

PARAMETER {
    gCl = 0.0013 (mho/cm2)            
    ECl = -20 (mV)                
    Cahalf_mCl = 0.0015 (mM)      
}

ASSIGNED {
    v (mV)                         
    cai (mM)                       
    icl (mA/cm2)                   
    mCl                            
}

BREAKPOINT {
    mCl = mCl_RPR(cai)             : calcium concentration
    icl = gCl * mCl * (v - ECl)   
}

INITIAL {
    mCl = 0.5
}

FUNCTION mCl_RPR(ca) {
    mCl_RPR = 1.0 / (1.0 + (Cahalf_mCl / ca) * (Cahalf_mCl / ca) * (Cahalf_mCl / ca) * (Cahalf_mCl / ca))
}
