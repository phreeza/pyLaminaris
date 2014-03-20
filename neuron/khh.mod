TITLE khh.mod Version 3 Hodgkin-Huxley squid K channel
 
COMMENT
 written by Jonathan Z Simon, jzsimon@isr.umd.edu
 but core code taken from Hines & Carnivale, "Expanding NEURON with NMODL"
 except adding more variables.
ENDCOMMENT

UNITS {
        (S)  = (siemens)
        (mV) = (millivolt)
        (mA) = (milliamp)
}
 
NEURON {
        SUFFIX khh
        USEION k READ ek WRITE ik
        RANGE  gkbar, gk, ik, alphaVHalf, alphaK, alpha0, betaVHalf, betaK, beta0, q10, T0
}

COMMENT
 Should be RANGE, not GLOBAL, or they cannot vary from cell to cell:
   alphaVHalf, alphaK, alpha0, betaVHalf, betaK, beta0, q10, T0
ENDCOMMENT
 
PARAMETER {
        gkbar = .036 (S/cm2)	<0,1e9>
        alpha0 = 0.1 (/ms)	    <0,1e9>
        alphaVHalf = -55 (mV)
        alphaK = 10.0 (mV)		<0,1e9>
        beta0 = 0.125 (/ms)	    <0,1e9>
        betaVHalf = -65 (mV)
        betaK = 80.0 (mV)	    <0,1e9>
        q10	= 3.0       	    <0,1e9>
        T0	= 6.3 (degC)
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
        gk = gkbar*n*n*n*n
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
         LOCAL x
         x = (Vm-alphaVHalf)*alphaKInv
         if (fabs(x) > 1e-6) {            : Traps for 0 in denominator
                alpha = alpha0*x/(1-exp(-x))
        }else{
                alpha = alpha0/(1 - 0.5*x)
        }
}
 
FUNCTION beta(Vm (mV)) (/ms) { 
         beta = beta0*exp((-Vm+betaVHalf)*betaKInv)
}
