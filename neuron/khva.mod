TITLE khva.mod   High Voltage Activation (~ Kv1.2? Kv3.1?) potassium channels
 
COMMENT
 written by Jonathan Z Simon, jzsimon@isr.umd.edu
 but core code taken from Hines & Carnivale, "Expanding NEURON with NMODL"
 except adding more variables
 Parameters & Functions taken from Rathouz & Trussel,
  J. Neurophysiol. 80: 2824-2835, 1998. (NM, Chick, Slice)
 klva.mod & khva.mod are identical except for suffix & defaults in PARAMETER
ENDCOMMENT
 
UNITS {
        (S)  = (siemens)
        (mV) = (millivolt)
        (mA) = (milliamp)
}
 
NEURON {
        SUFFIX khva
        USEION k READ ek WRITE ik
        RANGE  gkbar, gk, ik, alphaVHalf, alphaK, alpha0, betaVHalf, betaK, beta0, q10, T0
}

COMMENT
 Should be RANGE, not GLOBAL, or they cannot vary from cell to cell:
   alphaVHalf, alphaK, alpha0, betaVHalf, betaK, beta0, q10, T0
ENDCOMMENT
 
PARAMETER {
        gkbar = 0.04 (S/cm2)	<0,1e9>
        alpha0 = 0.11 (/ms)		<0,1e9>
        alphaVHalf = -19 (mV)
        alphaK = 9.10 (mV)		<0,1e9>
        beta0 = 0.103 (/ms)		<0,1e9>
        betaVHalf = -19 (mV)
        betaK = 20.0 (mV)		<0,1e9>
        q10	= 2.0				<0,1e9>
        T0	= 23 (degC)			<0,1e9>
}
 
ASSIGNED {
        gk (mho/cm2)
        ik (mA/cm2)
        ratefac
        alphaKInv (/mV)
        betaKInv (/mV)
        v (mV)
        celsius (degC)
        ek (mV)
}

STATE { n }
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gkbar*n : linear in n
        ik = gk*(v - ek)      
}
 
INITIAL {
        n = alpha(v)/(alpha(v) + beta(v)) : from states with n' vanishing
        ratefac = q10^((celsius - T0)/10(degC)) : can't vary with time
        alphaKInv = 1/alphaK : multiply by inverse is faster than divide
        betaKInv = 1/betaK
}

DERIVATIVE states {
        n' = ((1-n)*alpha(v) - n*beta(v))*ratefac
}
 
FUNCTION alpha(Vm (mV)) (/ms) { 
       alpha = alpha0*exp((Vm-alphaVHalf)*alphaKInv) 
}
 
FUNCTION beta(Vm (mV)) (/ms) { 
        beta  = beta0*exp((-Vm+betaVHalf)*betaKInv)
}
