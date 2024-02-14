TITLE CLSTmainen.mod   squid sodium, potassium, and chloride leak channels

COMMENT
This is a modification of hh.mod, the original Hodgkin-Huxley model provided by NEURON.
I have added Na-K pumps and the ability to left-shift a portion of the sodium channels.
These modifications reproduce the dynamics of Boucher PA, Joos B, Morris CE (2012) 'Coupled left-shift of Nav channels: modeling the Na(+)-loading and dysfunctional excitability of damaged axons' J Comput Neurosci,
making this model available for widespread use.

Further, I have included the NaV1.2, NaV1.6, and KV kinetics from Hu et al. (2009), which were originally from Mainen and Sejnowski (1996).

The following comments (in quotes) are those of the original hh.mod author:

"This is the original Hodgkin-Huxley treatment for the set of sodium, 
potassium, and leakage channels found in the squid giant axon membrane.
("A quantitative description of membrane current and its application 
conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
Membrane voltage is in absolute mV and has been reversed in polarity
from the original HH convention and shifted to reflect a resting potential
of -65 mV.
Remember to set celsius=6.3 (or whatever) in your HOC file.
See squid.hoc for an example of a simulation using this model.
SW Jaslove  6 March, 1992"
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
       SUFFIX CLSTmainen
       USEION na READ ena, nai, nao WRITE ina CHARGE 1
       USEION k READ ek, ki, ko WRITE ik CHARGE 1
       USEION cl READ ecl, cli, clo WRITE icl CHARGE -1

       :NONSPECIFIC_CURRENT il
       ::: vLS = Voltage left-shift VARIABLE (Cathy Morris and Bela Joos)
       RANGE ena0, ek0
       RANGE gkbar, gk, gkl, gl, el, INaKmax, ink, ink_last, Kmko, Kmnai, nai0, nao0, ki0, ko0, time0, ACpotassium, qPump_0, qNa_0, qK_0, qGate_0, tauRAMP, RAMP
       RANGE gnabar, gna, gnal, m, h, n, mLS, hLS, nLS, minf, hinf, ninf, mtau, htau, ntau, mLSinf, hLSinf, nLSinf, mLStau, hLStau, nLStau
       RANGE gnabar2, gna2, gnal2, m2, h2, mLS2, hLS2, minf2, hinf2, mtau2, htau2, mLSinf2, hLSinf2, mLStau2, hLStau2
       RANGE I_Nav, I_Nav2
       RANGE ChargeNav, ChargeNav2, ChargeNavTot
       RANGE gnabarTotal, Tref
       RANGE vLS, vLS0, vLeftShift, AC, timeLS
       RANGE vLS2, vLS02, vLeftShift2, AC2 :<---Properties of Nav1.2
       RANGE vRS, vRS0 :<---Properties of Nav1.2
       RANGE vShiftNa, vShiftK, vShift_h
       RANGE stimON, vThreshold, spike, gnabarSpike, tspike, tup, tdown, vpeak, tpeak
       RANGE slope :<---Dummy variable
       GLOBAL steepen
       GLOBAL Ra, Rb, Rd, Rg, tha, thi1, thi2, qa, qi, qinf
       GLOBAL qa2
       GLOBAL Ra_n, Rb_n, tha_n, qa_n
       GLOBAL tadj
       RANGE seg_len, seg_area, area_per_unit_length
       RANGE availability, activation
       RANGE taus, mtaus
       RANGE gnav
       RANGE disableRightShift_minf2, disableRightShift_mtau2
       RANGE disableRightShift_hinf2, disableRightShift_htau2
       RANGE vRS03, vRS3
       RANGE scaleNaV2
       
  THREADSAFE : assigned GLOBALs will be per thread
}

PARAMETER {
        disableRightShift_minf2 = 0
        disableRightShift_mtau2 = 0
        disableRightShift_hinf2 = 0
        disableRightShift_htau2 = 0
        scaleNaV2 = 1.0 
        

        Tref = 37.0
        vShiftNa = 0.0 (mV)  : Voltage shift applied to all SODIUM channels (calibration)
        vShiftK = 0.0 (mV)    : Voltage shift applied to all POTASSIUM channels (calibration)
        vShift_h = 0.0 (mV)

        timeLS = 200.0 (ms)

        vThreshold = 0.0 (mV)
        vpeak = -1000.0 (mV)
        gnabarSpike = 0.0 (S/cm2) <0,1e9>

        stimON = -1.0 (ms)
        tpeak = -1.0 (ms)
        tup = -1.0 (ms)
        tdown = -1.0 (ms)
        tspike = -1.0 (ms)
        spike = 0.0 
        slope = 0.0
        CLSpotassium = 0 (mV) :::Left-Shift potassium channels as well?

        time0 = 0.0 (ms)

        vLS0 = 0.0 (mV)
        vLeftShift = 0.0 (mV)  :Left-shift applied to damaged portion of Nav1.6 channels (AC)

        vLS02 = 0.0 (mV)
        vLeftShift2 = 0.0 (mV) :Left-shift applied to damaged portion of Nav1.2 channels (AC2)

        steepen = 1.00

        vRS0 = 0.0 (mV) : 13.0 (mV) #-6.5  :::Relative right-shift of Nav1.2
        vRS03 = 0.0 (mV)

        AC = 0.0  <0,1.0> ::: Proportion of affected (left-shifted) NaV1.6 channels on node
        AC2 = 0.0  <0,1.0> ::: Proportion of affected (left-shifted) NaV1.2 channels on node
        ACpotassium = 0.0  <0,1.0> ::: Proportion of affected (left-shifted) POTASSIUM channels on node
        

        INaKmax = 9.09e-2       (mA/cm2) <0,1e6>



        ena0 =  51.5 (mV)
        ek0 =  -81.3 (mV)
        el = -59.9 (mV)

        ink0 = 0.010731         (mA/cm2) <0,1e6> 
        Kmnai =    10          (mM)    <0,1e6>
        Kmko =     3.5         (mM)    <0,1e6>



        :::PAB
        nai0 = 20.0 (mM)
        nao0 = 154.0 (mM) 
        ki0 = 150.0 (mM)
        ko0 = 6.0 (mM)

        cli0 = 11.65 (mM) :11.5*mM :9.9*mM
        clo0 = 125.0 (mM) :124.5*mM :123.27*mM


        gnabar2 = 0.0 (S/cm2) <0,1e9> :<--- Nav1.2
        gnabar = 0.12 (S/cm2) <0,1e9>
        gkbar = 0.036 (S/cm2) <0,1e9>
        gl = 0.0005 (S/cm2) <0,1e9>
        gnal = 0.00025 (S/cm2) <0,1e9> :0.00025 <--PAB
        gkl = 0.0001 (S/cm2) <0,1e9> :0.0001 <--PAB


        

        :::Temperature Factors ("Q-ten's" as they say)
        qPump_0 = 1.9
        qNa_0 = 1.4
        qK_0 = 1.1
        qGate_0 = 3.0

        tauRAMP = 3000000

        seg_len = 0.0 (um)



        :Mainen1996 param names
        tha  = -43  (mV)        : v 1/2 for act     
        qa   = 6    (mV)        : act slope 
        qa2   = 7    (mV)        : act slope     
        Ra   = 0.182    (/ms)       : open (v)      
        Rb   = 0.124    (/ms)       : close (v)     

        thi1  = -50 (mV)        : v 1/2 for inact   
        thi2  = -75 (mV)        : v 1/2 for inact   
        qi   = 5    (mV)            : inact tau slope
        thinf  = -72    (mV)        : inact inf slope   
        qinf  = 6.2 (mV)        : inact inf slope
        Rg   = 0.0091   (/ms)       : inact (v) 
        Rd   = 0.024    (/ms)       : inact recov (v)

        tha_n  = 25  (mV)    : v 1/2 for act   (-42)
        qa_n   = 9  (mV)    : act slope   
        Ra_n   = 0.02 (/ms)   : max act rate   
        Rb_n   = 0.002  (/ms)   : max deact rate 
 

}

STATE {
      time (ms)

      m
      h
      n
      mLS
      hLS
      nLS

      m2
      h2
      mLS2
      hLS2

      ChargeNav
      ChargeNav2
      ChargeNavTot

      

      
}

ASSIGNED {        
        gnabarTotal (S/cm2) :<---Combined gnabar+gnabar2 (mostly to avoid ZeroDivisionErrors in data analysis)
        v (mV)
        vLS (mV)

        vRS (mV)
        vRS3 (mV)
        vLS2 (mV)

        celsius (degC)
        ena (mV)
        ek (mV)
        ecl (mV)

        nai (mM)
        nao (mM)
        ki (mM)
        ko (mM)
        cli (mM)
        clo (mM)

        ink (mA/cm2)
        ink_last (mA/cm2)

        gna (S/cm2)
        gk (S/cm2)
        ina (mA/cm2)
        I_Nav (mA/cm2)
        I_Nav2 (mA/cm2)

        seg_area (um2)
        area_per_unit_length (um2)

        ik (mA/cm2)
        icl (mA/cm2)
        minf
        hinf
        ninf
        mtau (ms) 
        htau (ms) 
        ntau (ms)

        mLSinf
        hLSinf
        nLSinf
        mLStau (ms) 
        hLStau (ms) 
        nLStau (ms) 

        taus (ms)
        mtaus (ms)
        availability
        activation
        gnav (S/cm2)

        gna2 (S/cm2)
        minf2
        hinf2
        mtau2 (ms) 
        htau2 (ms) 

        mLSinf2
        hLSinf2
        mLStau2 (ms) 
        hLStau2 (ms) 


        diam (um)
        L (um)
        RAMP
        area (um2)
}

INITIAL {
  v = -70.0 :-59.8137886
  gnabarTotal = gnabar + gnabar2
  vLS = vLS0
  vLS2 = vLS0
  vRS = vRS0
  vRS3 = vRS03
  time = time0 


  ratesMainen96(v)
  m = minf
  h = hinf
  n = ninf
  mLS = mLSinf
  hLS = hLSinf
  nLS = nLSinf


  m2 = minf2
  h2 = hinf2
  mLS2 = mLSinf2
  hLS2 = hLSinf2

  ChargeNav = 0.0
  ChargeNav2 = 0.0
  ChargeNavTot = 0.0

  availability = (gnabar*hinf + gnabar2*hinf2)/(gnabar + gnabar2)
  activation = (gnabar*minf + gnabar2*minf2)/(gnabar + gnabar2)


  RAMP = 0.0
}

: currents
BREAKPOINT {
  LOCAL qPump, qNa, qK

 RAMP = qPump_0
 seg_area = seg_len*3.1415926536*diam
 area_per_unit_length = 3.1415926536*diam



  qPump = qPump_0^((celsius - Tref)/10.0)
  qNa = qNa_0^((celsius - Tref)/10.0)
  qK = qK_0^((celsius - Tref)/10.0)
  SOLVE states METHOD cnexp 
  :ink_last = ink
  ink= qPump*INaKmax/(((1 + (Kmnai/nai))^3)*((1 + Kmko/ko)^2))

  gna = qNa*gnabar*(m*m*m*h*(1.0 - AC) + mLS*mLS*mLS*hLS*AC)
  gna2 = scaleNaV2*qNa*gnabar2*(m2*m2*m2*h2*(1.0 - AC2) + mLS2*mLS2*mLS2*hLS2*AC2)


  gk = qK*gkbar*(n*(1.0 - ACpotassium) + nLS*ACpotassium )  :<---Mainen


  I_Nav = gna*(v - ena)*area_per_unit_length :*area_per_unit_length  :*seg_area  :*area*(1e-8)
  I_Nav2 = gna2*(v - ena)*area_per_unit_length :*area_per_unit_length :*seg_area  :*area*(1e-8)

  ina = (gna + gna2 + gnal)*(v - ena) + 3*ink
  
  ik = (gk + gkl)*(v - ek) - 2*ink 
  icl = gl*(v - ecl)

  availability = (gnabar*h + gnabar2*h2)/(gnabar + gnabar2)
  activation = (gnabar*m + gnabar2*m2)/(gnabar + gnabar2)
  gnav = (gnabar*gna + gnabar2*gna2)/(gnabar + gnabar2)
}



: states
DERIVATIVE states {
        LOCAL nai_prime, ki_prime, qGate
        qGate = qGate_0^((celsius - Tref)/10.0)

        time' = 1.0

        :rates(v)
        ratesMainen96(v)
        m' =  qGate*(minf-m)/mtau
        h' = qGate*(hinf-h)/htau
        n' = qGate*(ninf-n)/ntau

        mLS' =  qGate*(mLSinf-mLS)/mLStau
        hLS' = qGate*(hLSinf-hLS)/hLStau
        nLS' =  qGate*(nLSinf-nLS)/nLStau

        m2' =  qGate*(minf2-m2)/mtau2
        h2' = qGate*(hinf2-h2)/htau2

        mLS2' =  qGate*(mLSinf2-mLS2)/mLStau2
        hLS2' = qGate*(hLSinf2-hLS2)/hLStau2


      ChargeNav' = I_Nav
      ChargeNav2' = I_Nav2
      ChargeNavTot' =  I_Nav+I_Nav2


}


UNITSOFF
PROCEDURE ratesMainen96(vm (mV)) {  
    LOCAL temp, q10, alpha, beta, alpha2, beta2, vshift, vshift2, vshift3


    vshift = 8.0 :(mV)
    vshift2 = 8.0 - vRS0 :(mV)
    vshift3 = 8.0 - vRS03 :(mV)
    temp = 23   (degC)      : original temp 
    q10  = 2.3          : temperature sensitivity
    tadj = q10^((celsius - temp)/10)
    :tadj = 1.0

    if (t < stimON) {
    ChargeNav = 0.0
    ChargeNav2 = 0.0
    ChargeNavTot = 0.0  
    }



    alpha = AlphaV(vm+vshift,tha,Ra,qa)
    beta = BetaV(vm+vshift,tha,Rb,qa)
    mtau = 1/tadj/(alpha+beta)
    minf = alpha/(alpha+beta)

    :"h" inactivation 
    alpha = AlphaV(vm+vshift,thi1,Rd,qi)
    beta = BetaV(vm+vshift,thi2,Rg,qi)
    htau = 1/tadj/(alpha+beta)
    hinf = 1/(1+exp((vm+vshift-thinf)/qinf))

    :"n" Potassium activation
    alpha = AlphaV(vm,tha_n,Ra_n,qa_n)
    beta = BetaV(vm,tha_n,Rb_n,qa_n)
    ntau = 1/tadj/(alpha+beta)
    ninf = alpha/(alpha+beta)


    :"m2" second sodium activation system
    if (disableRightShift_minf2 == 1) {
            alpha2 = AlphaV(vm+vshift3,tha,Ra,qa2)
            beta2 = BetaV(vm+vshift3,tha,Rb,qa2)
            minf2 = alpha2/(alpha2+beta2)

    }else{
            alpha2 = AlphaV(vm+vshift2,tha,Ra,qa2)
            beta2 = BetaV(vm+vshift2,tha,Rb,qa2)
            minf2 = alpha2/(alpha2+beta2)
    }


    if (disableRightShift_mtau2 == 1) {
            alpha2 = AlphaV(vm+vshift3,tha,Ra,qa2)
            beta2 = BetaV(vm+vshift3,tha,Rb,qa2)
            mtau2 = 1/tadj/(alpha2+beta2)

    }else{
            alpha2 = AlphaV(vm+vshift2,tha,Ra,qa2)
            beta2 = BetaV(vm+vshift2,tha,Rb,qa2)
            mtau2 = 1/tadj/(alpha2+beta2)
    }

    :"h2" second sodium inactivation system
    if ( disableRightShift_hinf2 == 1) {
            hinf2 = 1/(1+exp((vm+vshift3-thinf)/qinf))
    }else{
            hinf2 = 1/(1+exp((vm+vshift2-thinf)/qinf))
    }

    if ( disableRightShift_htau2 == 1) {

            alpha2 = AlphaV(vm+vshift3,thi1,Rd,qi)
            beta2 = BetaV(vm+vshift3,thi2,Rg,qi)
            htau2 = 1/tadj/(alpha2+beta2)
    }else{
            alpha2 = AlphaV(vm+vshift2,thi1,Rd,qi)
            beta2 = BetaV(vm+vshift2,thi2,Rg,qi)
            htau2 = 1/tadj/(alpha2+beta2)

            
    }


    if ( t <= timeLS) {
               vLS = vLeftShift*(   0.5*( 1 + tanh((7.0/timeLS)*t - 3.0) )   )
               vLS2 = vLeftShift2*(   0.5*( 1 + tanh((7.0/timeLS)*t - 3.0) )   )
       }
       if (t > timeLS) {
           vLS = vLeftShift
           vLS2 = vLeftShift2

       }

        if (t > stimON) {
               if ( v > vThreshold ) {
                   if ( spike==0.0 ) { 
                       tup = t
                       gnabarSpike = gnabar + gnabar2
                       } 
                   spike = 1.0 
                   }

              if (tup != -1.0) {
                   if (t > tup) {
                      if (v < vThreshold) {
                          if (tdown==-1.0) {
                              tdown = t
                              
                              }
                          }
                          tspike = (tdown+tup)/2.0
                       }
                  }
              :if (spike == 1.0) {
                  if (tdown==-1.0) {
                      if (v > vpeak) {
                          vpeak = v
                          tpeak = t
                          }
                      }
              :     }
              }


    :"mLS" Affected (Left Shifted) sodium activation system
    alpha = AlphaV(vm+vshift+vLS,tha,Ra,qa)
    beta = BetaV(vm+vshift+vLS,tha,Rb,qa)
    mLStau = 1/tadj/(alpha+beta)
    mLSinf = alpha/(alpha+beta)

    :"hLS" Affected (Left Shifted) sodium inactivation system
    alpha = AlphaV(vm+vshift+vLS,thi1,Rd,qi)
    beta = BetaV(vm+vshift+vLS,thi2,Rg,qi)
    hLStau = 1/tadj/(alpha+beta)
    hLSinf = 1/(1+exp((vm+vshift+vLS-thinf)/qinf))

    :"nLS" Affected (Left Shifted) potassium activation system
    alpha = AlphaV(vm+vLS,tha_n,Ra_n,qa_n)
    beta = BetaV(vm+vLS,tha_n,Rb_n,qa_n)
    nLStau = 1/tadj/(alpha+beta)
    nLSinf = alpha/(alpha+beta)


    :"mLS2" second Affected (Left Shifted) sodium activation system
    alpha2 = AlphaV(vm+vshift2+vLS2,tha,Ra,qa2)
    beta2 = BetaV(vm+vshift2+vLS2,tha,Rb,qa2)
    mLStau2 = 1/tadj/(alpha2+beta2)
    mLSinf2 = alpha2/(alpha2+beta2)

    :"hLS2" second Affected (Left Shifted) sodium inactivation system
    alpha2 = AlphaV(vm+vshift2+vLS2,thi1,Rd,qi)
    beta2 = BetaV(vm+vshift2+vLS2,thi2,Rg,qi)
    hLStau2 = 1/tadj/(alpha2+beta2)
    hLSinf2 = 1/(1+exp((vm+vshift2-thinf)/qinf))

    taus = (gnabar*htau + gnabar2*htau2)/(gnabar+gnabar2)
    mtaus = (gnabar*mtau + gnabar2*mtau2)/(gnabar+gnabar2)




}

FUNCTION AlphaV(s,th,A,q) {
  AlphaV = A * q * efun((th - s)/q)
} 

FUNCTION BetaV(s,th,A,q) {
  BetaV = A * q * efun((s - th)/q)
} 
UNITSON


FUNCTION efun(z) {
  if (fabs(z) < 1e-6) {
    efun = 1 - z/2
  }else{
    efun = z/(exp(z) - 1)
  }
}

