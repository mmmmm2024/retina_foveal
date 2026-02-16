COMMENT
  IinjLT.mod
  Generates a train of photocurrent injections
  User specifies duration of each Iinj, interpulse interval (ton and toff),
  and number of pulses.
  6/30/2003 RARE LAB
           ___            ___            _ _
          /   \          /   \            |
 delay    |    \   toff  |    \           | amp
__________|     \________|     \_______  _|_
            ton
            
   num:      1              2          ...
ENDCOMMENT


NEURON {
	POINT_PROCESS IinjLT_oncone
	RANGE del, ton, toff, num, amp,ssI,i
	ELECTRODE_CURRENT i
}

UNITS {
	(pA) = (picoamp)
        (nA) = (nanoamp)
}

PARAMETER {
	del  = 500 (ms)
	ton  = 2000 (ms) <0, 1e9>	: duration of ON phase
	toff = 6000 (ms) <0, 1e9>	: duration of OFF phase
	num  = 1			: how many to deliver
	amp  = 0 (pA)	  : absolute amplitude of the light induced current
      ssI  = 40 (pA)     : steady-state current (dark current)
        
}

ASSIGNED {
        Ncount     : counter of the number of the flashes/injections
	ival (nA)
	i (nA)
	on
	tally			: how many more to deliver
	tr (ms)   : the relative time in each flash 
        Part1
        Part2
        Part3
        ssInA (nA)
        ampnA (nA)
}

INITIAL {
	i = 0
        ssInA=0.001*ssI
        ampnA=amp*0.001
	ival = 0
	tally = num
        Ncount=0
	if (tally > 0) {
		net_send(del, 1)
		on = 0
		tally = tally - 1
	}
}

:only 1(0 ~ 5000[ms]) trial
BREAKPOINT{
	if(t<=2000){
		i = 0
	}

	if(t>2000 && t<3000){
	i = 0.0028*log(2*t-3999) + 0
	}

	if(t>=3000){
		i = 1 / (1.5*t-4450) + 0
	}
}

NET_RECEIVE (w) {
	: ignore any but self-events with flag == 1
	if (flag == 1) {
		if (on == 0) {
			: turn it on
                        Ncount=Ncount+1
			ival = ampnA
			on = 1
			: prepare to turn it off
			net_send(ton, 1)
		} else {
			: turn it off
			ival = 0
			on = 0
			if (tally > 0) {
				: prepare to turn it on again
				net_send(toff, 1)
				tally = tally - 1
			}
		}
	}
}
