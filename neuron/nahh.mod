TITLE nahh.mod Version 3 Hodgkin-Huxley squid Na channel
 
COMMENT
 written by Jonathan Z Simon, jzsimon@isr.umd.edu
 but core style taken from Hines & Carnivale, "Expanding NEURON with NMODL"
 except adding more variables.
ENDCOMMENT
 
UNITS {
        (S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX nahh
        USEION na READ ena WRITE ina
        RANGE  gnabar, gna, ina, alphamVHalf, alphamK, alpham0, betamVHalf, betamK, betam0, alphahVHalf, alphahK, alphah0, betahVHalf, betahK, betah0, q10, T0
}

COMMENT
 Should be RANGE, not GLOBAL, or they cannot vary from cell to cell:
   alphaVHalf, alphaK, alpha0, betaVHalf, betaK, beta0, q10, T0
ENDCOMMENT
 
PARAMETER {
        gnabar = 0.12 (S/cm2)	<0,1e9>
        alpham0 = 1.0 (/ms)		<0,1e9>
        alphamVHalf = -40 (mV)
        alphamK = 10.0 (mV)		<0,1e9>
        betam0 = 4.0 (/ms)		<0,1e9>
        betamVHalf = -65 (mV)
        betamK = 18.0 (mV)		<0,1e9>
        alphah0 = 0.07 (/ms)	<0,1e9>
        alphahVHalf = -65 (mV)
        alphahK = 20.0 (mV)		<0,1e9>
        betah0 = 1.0 (/ms)		<0,1e9>
        betahVHalf = -35 (mV)
        betahK = 10.0 (mV)		<0,1e9>
        q10	= 3.0				<0,1e9>
        T0	= 6.3 (degC)		<0,1e9>
}
 
STATE { m h }
 
ASSIGNED {
		gna (mho/cm2)
		ina (mA/cm2)
		ratefac
        alphamKInv (/mV)
		betamKInv (/mV)
        alphahKInv (/mV)
		betahKInv (/mV)
        v (mV)
        celsius  (degC)
        ena (mV)
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
        ina = gna*(v - ena)      
}
 
INITIAL {
        m = alpham(v)/(alpham(v) + betam(v)) : from states with m' vanishing
        h = alphah(v)/(alphah(v) + betah(v)) : from states with h' vanishing
		ratefac = q10^((celsius - T0)/10(degC)) : can't vary with time
		alphamKInv = 1/alphamK : multiply by inverse is faster than divide
		betamKInv = 1/betamK
		alphahKInv = 1/alphahK
		betahKInv = 1/betahK
}

 
DERIVATIVE states { 
		m' = ((1-m)*alpham(v) - m*betam(v))*ratefac
		h' = ((1-h)*alphah(v) - h*betah(v))*ratefac
}
 
FUNCTION alpham(Vm (mV)) (/ms) { 
		 LOCAL x
		 x = (Vm-alphamVHalf)*alphamKInv
		 if (fabs(x) > 1e-6) {            : Traps for 0 in denominator
                alpham = alpham0*x/(1-exp(-x))
        }else{
                alpham = alpham0/(1 - 0.5*x)
        }
}
 
FUNCTION betam(Vm (mV)) (/ms) { 
         betam = betam0*exp((-Vm+betamVHalf)*betamKInv)
}
 
FUNCTION alphah(Vm (mV)) (/ms) { 
         alphah = alphah0*exp((-Vm+alphahVHalf)*alphahKInv)
}
 
FUNCTION betah(Vm (mV)) (/ms) { 
         betah = betah0/(exp((-Vm+betahVHalf)*betahKInv) + 1)
}
 
