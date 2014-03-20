TITLE Sinusoid Conductance
                 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  SUFFIX sing
  POINT_PROCESS SinConductance
  RANGE onset, gAC,gDC, e, i, g
  NONSPECIFIC_CURRENT i
}

UNITS {
  (nA)  = (nanoamp)
  (mV)  = (millivolt)
  (nS)  = (nanomho)
}

PARAMETER {
  onset = 5  (ms)
  T  = 0.25     (ms)
  gDC = 11.88 (nS)
  gAC = 6.0 (nS)
  e= 0.0    (mV)
  v   (mV)
}

ASSIGNED { 
  i     (nA)  
  g     (nS)
}

UNITSOFF

BREAKPOINT { LOCAL tt
  tt= (t-onset)
  if (t>onset && tt<1630) { 
    g = gAC * sin(tt/T*2*3.14) + gDC
  }
  else {g = 0}
  i = g * (v - e)
}
UNITSON
