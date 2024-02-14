:Reference :Colbert and Pan 2002

NEURON	{
	SUFFIX NaTa6_t
	USEION na READ ena WRITE ina
	RANGE vRS0, vRS03, gNaTa6_tbar, gNaTa6_t, ina
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gNaTa6_tbar = 0.00001 (S/cm2)
	vRS0 = -6.0 (mV) : 0.0 (mV)
	vRS03 = 0.0 (mV)
}

ASSIGNED	{
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTa6_t	(S/cm2)
	mInf
	mTau    (ms)
	mAlpha
	mBeta
	hInf
	hTau    (ms)
	hAlpha
	hBeta
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gNaTa6_t = gNaTa6_tbar*m*m*m*h
	ina = gNaTa6_t*(v-ena)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((34-21)/10)
	
  UNITSOFF
    if(v == -(38-vRS0)){
    	v = v+0.0001
    }
		mAlpha = (0.182 * (v- -(38-vRS0)))/(1-(exp(-(v- -(38-vRS0))/6)))
		mBeta  = (0.124 * (-v -(38-vRS0)))/(1-(exp(-(-v -(38-vRS0))/6)))
		mTau = (1/(mAlpha + mBeta))/qt
		mInf = mAlpha/(mAlpha + mBeta)

    if(v == -(66-vRS0)){
      v = v + 0.0001
    }

		hAlpha = (-0.015 * (v- -(66-vRS0)))/(1-(exp((v- -(66-vRS0))/6)))
		hBeta  = (-0.015 * (-v -(66-vRS0)))/(1-(exp((-v -(66-vRS0))/6)))
		hTau = (1/(hAlpha + hBeta))/qt
		hInf = hAlpha/(hAlpha + hBeta)
	UNITSON
}
