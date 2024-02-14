:Reference :Colbert and Pan 2002
:comment: took the NaTa and shifted both activation/inactivation by 6 mv

NEURON	{
	SUFFIX NaTa2_t
	USEION na READ ena WRITE ina
	RANGE vRS0, vRS03, DeltaVrs, gNaTa2_tbar, gNaTa2_t, ina
	RANGE scaleNavTwo
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gNaTa2_tbar = 0.00001 (S/cm2)
	vRS0 = 6.0 (mV)
	vRS03 = 0.0 (mV)
	DeltaVrs = 0.0 (mV)
	scaleNavTwo = 1.0
}

ASSIGNED	{
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTa2_t	(S/cm2)
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
	gNaTa2_t = scaleNavTwo*gNaTa2_tbar*m*m*m*h
	ina = gNaTa2_t*(v-ena)
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
    if(v == -(38.0 - (vRS0+DeltaVrs))){
    	v = v+0.0001
    }
		mAlpha = (0.182 * (v- -(38.0 - (vRS0+DeltaVrs))))/(1-(exp(-(v- -(38.0 - (vRS0+DeltaVrs)))/6)))
		mBeta  = (0.124 * (-v -(38.0 - (vRS0+DeltaVrs))))/(1-(exp(-(-v -(38.0 - (vRS0+DeltaVrs)))/6)))
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = (1/(mAlpha + mBeta))/qt

    if(v == -(66.0 - (vRS0+DeltaVrs))){
      v = v + 0.0001
    }
		hAlpha = (-0.015 * (v- -(66.0 - (vRS0+DeltaVrs))))/(1-(exp((v- -(66.0 - (vRS0+DeltaVrs)))/6)))
		hBeta  = (-0.015 * (-v -(66.0 - (vRS0+DeltaVrs))))/(1-(exp((-v -(66.0 - (vRS0+DeltaVrs)))/6)))
		hInf = hAlpha/(hAlpha + hBeta)
		hTau = (1/(hAlpha + hBeta))/qt
	UNITSON
}
