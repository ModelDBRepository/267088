TITLE peak.mod

COMMENT
This mechanism stores the peak voltage at every segment
ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mol) = (1)
    (molar) = (mol/liter)
    (mM) = (millimolar)
    (uM) = (micromolar)
    (um) = (micrometer)
    FARADAY = (faraday) (coulombs)
}

NEURON {
    SUFFIX peak
    RANGE vpeak
    RANGE tpeak
}

PARAMETER {
    vpeak = -1000.0 (mV)
    tpeak = -1000.0 (ms)
}

BREAKPOINT {
    if (v > vpeak) {
      vpeak = v
      tpeak = t
      }
                      
}