:yonemoto simPR.c
NEURON {
    POINT_PROCESS RPRInputoff
    :SUFFIX RPRInput
    RANGE del, amp, i, v, ton, toff , N, num      
    ELECTRODE_CURRENT i    
}

UNITS {
	(pA) = (picoamp)
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uF) = (microfarad)
}

PARAMETER {
    :i (nA)
    del = 2000 (ms)
    ton =1000 (ms) <0, 1e9>
    toff =0 (ms) <0, 1e9>              
    amp = 0 (pA)
    area = 200 (um2)
    num  = 1
    :cm = 10 (uF)                  
}

ASSIGNED {   
    Ncount                      
    i (nA)
    ival (nA)
    :onset (ms)                        
    :t (ms)
    on
    ampnA (nA)
    :light  (ms)
    tal
    light (ms)                       
}

INITIAL {
    i = 0             : Initial
    ampnA = amp * 0.001
    ival = 0
    	tal = num
        Ncount =0
	if (tal > 0) {
        on = 1
		net_send(del, 1)
		tal = tal - 1
	}
}

BREAKPOINT {
    LOCAL photo
    
    light = t - del - (ton + toff) * (Ncount-1)
    photo = (32*( 1-exp(- (light/1000 )/0.05  ) ) -33/(   1+exp(-   (  (light/1000) -3.8   )/0.45    ) ) +1-exp(  - (light/1000)  /0.8 ))/33.0
    if (on==1) {
        i = 0.04 - ival * photo 
    } else {
        i = 0.04 + ival
    }

}

NET_RECEIVE (w) {
	: ignore any but self-events with flag == 1
	if (flag == 1) {
		if (on == 0) {
			: turn it on
            ival = ampnA
			on = 1
			: prepare to turn it on
			if (tal > 0) {
				: prepare to turn it on again
				net_send(ton, 1)
				tal = tal - 1
			}
		} else {
			: turn it off
            Ncount=Ncount+1
            ival = 0
			on = 0
			net_send(toff, 1)
		}
	}
}