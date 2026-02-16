NEURON {
    POINT_PROCESS Iphoto_rod
    RANGE Idark, A, t1, t2, t3, b, amp, i
    ELECTRODE_CURRENT i
}

UNITS {
    (pA) = (picoamp)
		(nA) = (nanoamp)
}

PARAMETER {
    :t (ms):time 
    Idark = -40 (pA)          : Dark current
    t1 = 50 (ms)              : first term
    t2 = 450 (ms)             : second term
    t3 = 800 (ms)             : third term
    b = 3800 (ms)             
    amp = 0 (pA)             : photocurrent
}

ASSIGNED {
    i (nA)                    : Output
    IdarknA (nA)
    ampnA (nA)

}

INITIAL {
    i = 0
    IdarknA = Idark * 0.001             : Initial
    ampnA = amp * 0.001

}

BREAKPOINT {
    i = IdarknA + 0.01 * Photocurrent()
}

FUNCTION Photocurrent() {
    LOCAL Part1, Part2, Part3
    Part1 = (1 - exp(-t / t1))
    Part2 = 1 / (1 + exp(-(t - b) / t2))
    Part3 = (1 - exp(-t / t3))
    
    Photocurrent = Part1 - Part2 + Part3
}