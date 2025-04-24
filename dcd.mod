UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uF) = (microfarad)
	(nC) = (nanocoul)
	(F) = (farad)
	(Pa) = (N/meter2)
	(um) = (micron)
	(mol) = (mole)
	PI = (pi) (1)
}

NEURON {
  POINT_PROCESS DcDt
  THREADSAFE
  RANGE cm0, tbegin, tdur, dc, a, LJ_C, LJ_alpha, Delta, m, n, f, A, z0
  RANGE q, stm, rel_Zmin
  POINTER c
  NONSPECIFIC_CURRENT i
}
 
PARAMETER {
  cm0 = 1 (uF/cm2)
  tbegin = 10 (ms) 
  tdur = 1000 (ms)
  a = 32e-3 (um)
  dc (uF/cm2-ms)
  w = 0 (/ms)
  A = 1 (Pa)
  f = 1000 (/ms)
  zeR = 0 (um)
  oneR = 1 (um)
  z0 = 0 (um)
  u0 = 0 (um/ms)
  
  T = 309.15 (K)
  Rg = 8.314 (Pa-m3/mol-K)
  delta0 = 2.0e-3 (um)
  Delta = 1.2736e-3 (um)
  Delta_ = 1.4e-3 (um)
  pDelta = 1.0e5   (Pa)
  m = 3.9176    (1)
  n = 1.0134   (1)
  LJ_C = 18404.94148 (Pa)
  LJ_alpha = 1.6606e-3 (um)
  
  rhoL = 1075.0 (kg/m3)
  muL = 7.0e-4 (Pa-s)
  muS = 0.035 (Pa-s)
  kA = 0.24  (N/m)
  alpha = 7.56  (Pa-s)
  C0 = 0.62  (mol/m3)
  kH = 1.613e5   (Pa-m3/mol)
  P0 = 1.0e5  (Pa)
  Dgl = 3.68e-9    (m2/s)
  xi = 0.5e-3 (um)
  Zqs = 0.001e-3 (um)
  
  epsilon0 = 8.854e-12 (F/m)
  epsilonR = 1.0 (1)
  
  rel_Zmin = -0.49 (1)
}
 
ASSIGNED {
  c (uF/cm2)
  i (nA)
  v (mV)
  area (um2)
  Zmin (um)


  q (nC/cm2)
  stm (Pa)
}

STATE {
	ng (mol)
	U (um/ms)
	Z (um)
}


FUNCTION arealStrain(Z(um)) (1) {
	arealStrain = (Z/a) * (Z/a)
}

FUNCTION surface(Z(um)) (um2) {
	surface = PI * (a*a + Z*Z)
}

FUNCTION curvrad(Z(um)) (um) {
	if (Z==zeR) {
		curvrad = 1000 * a
	} else {
	curvrad = (a*a + Z*Z)/(2*Z)
	}
}

FUNCTION sign(i(um)) (1) {
	if (i<zeR) {
		sign = -1
	}
	if (i==zeR) {
		sign = 0
	} 
	if (i>zeR) {
		sign = 1
	}
}

FUNCTION PMavgPred (Z(um)) (Pa) {
	LOCAL LJ_factor
	if (Z<Zqs) {
		PMavgPred = 0
	} else {
		LJ_factor = LJ_alpha / (2 * Z + Delta)
		PMavgPred = LJ_C * (LJ_factor ^ m - LJ_factor ^ n)
	}
}

FUNCTION PEtot (Z(um)) (Pa) {
	PEtot = -(1e+6)*kA * arealStrain(Z) / curvrad(Z)
}

FUNCTION Pelec (Z(um), Qm(nC/cm2)) (Pa) {
	LOCAL rels, abs_perm, S0
	S0 = PI * a * a
	rels = S0 / surface(Z)
	abs_perm = epsilon0 * epsilonR
	Pelec = -(1e-10) * rels * Qm * Qm / (2 * abs_perm)
}

FUNCTION Pvtot (U(um/ms), R(um)) (Pa) {
	LOCAL PVleaflet, PVfluid
	PVleaflet = -12 * U * delta0 * muS / R^2
	PVfluid = -4 * U * muL / abs(R)
	Pvtot = (1000)*(PVleaflet + PVfluid)
}

FUNCTION volume (Z(um)) (um3) {
	volume = PI * a * a * Delta * (1 + (Z/(3*Delta)*(3+Z*Z/(a*a))))
}

FUNCTION gasmol2Pa (ng(mol), Vol(um3)) (Pa) {
	gasmol2Pa = ng * Rg * T / ((1e-18)*Vol)
}

FUNCTION gasmol2Paqs (Vol(um3)) (Pa) {
	gasmol2Paqs = P0 * PI * a * a * Delta / Vol
}

FUNCTION gasPa2mol (P(Pa), Vol(um3)) (mol) {
	gasPa2mol = P * (1e-18)*Vol / (Rg * T)
}

FUNCTION gasFlux (Z(um), Pg(Pa)) (mol/s) {
	gasFlux = (1e-6) * 2 * surface(Z) * Dgl * (C0-Pg/kH) / xi
}

FUNCTION abs(r(um)) (um) {
	if (r < 0) {
		abs = -r
	} else {
		abs = r
	}
}

FUNCTION sele(Z(um), ng(mol)) (mol/s) {
	if (Z<=Zqs) {
		sele = gasFlux(Z, gasmol2Paqs(volume(Z)))
	} else {
		sele = gasFlux(Z, gasmol2Pa(ng, volume(Z)))
	}
}

FUNCTION n_num(r(um2)) (1) {
	n_num = r / oneR^2
}
 
FUNCTION cm(t(ms), Z(um)) (uF/cm2) {
	LOCAL z2
	if (t > tbegin && t < (tbegin + tdur)) {
		
		if (Z==zeR) {
			cm = cm0
		} else {
			z2 = (a*a-Z*Z-Z*Delta)/(2*Z)
			cm = cm0*Delta/(a*a)*(Z+z2*log((2*Z+Delta)/Delta))
		} 
	} else {
		cm = cm0
	}
}
        
FUNCTION dcmdt(Z(um), U(um/ms))(uF/cm2-ms) {
	LOCAL ratio1, ratio2
	if (t > tbegin && t < (tbegin + tdur)) {
		ratio1 = (Z*Z + a*a) / (Z*(2*Z+Delta))
		ratio2 = (Z*Z + a*a) / (2 * Z*Z) * log((2*Z+Delta)/Delta)
		dcmdt = cm0*Delta/(a*a)*(ratio1-ratio2)*U
	} else {
		dcmdt = 0
	} 
}

FUNCTION dUdt(Z(um), U(um/ms), ng(mol)) (um/ms2) {
	LOCAL R, S, V, Pg, Pm, Pv, Pac, Ptot, accP, accNL
	
	if (t > tbegin && t < (tbegin + tdur)) {
	
		if (Z > a) {
			Z = a
		}
		if (Z < Zmin) {
			Z = Zmin
		}
		
		R = curvrad(Z)
		S = surface(Z)
		V = volume(Z)
		
		if (Z<=Zqs) {
			Pg = gasmol2Paqs(V)
		} else {
			Pg = gasmol2Pa(ng, V)
		}
		Pm = PMavgPred(Z)
		Pv = Pvtot(U, R)
		Pac = A * sin(2*PI*f*t)
		stm = Pac
		Ptot = Pm+Pg-P0-Pac+PEtot(Z)+Pv+Pelec(Z, c*v)
		accP = Ptot / (rhoL * R)
		accNL = -(3 * U ^ 2) / (2 * R)
		dUdt = accNL + (1e+6) * accP
	} else {
		dUdt = 0 (um/ms2)
		stm = 0
	}
}

FUNCTION Zbound(Z(um), t(ms)) (um) {
	if (t > tbegin && t < (tbegin + tdur)) {
		if (Z > a) {
			Z = a
		}
		if (Z < Zmin) {
			Z = Zmin
		}
	} else {
		Z = z0
	}
}

INITIAL {
	c = cm0
	Zmin = rel_Zmin * Delta
	ng = gasPa2mol(P0, PI * Delta * a* a)
	dc = 0
	U = u0
	Z = z0
	net_send(tbegin, 1)
}

BEFORE BREAKPOINT {
	Z = Zbound(Z, t)
	if (t<= tbegin || t >= (tbegin + tdur)) {
		ng = gasPa2mol(P0, PI * Delta * a* a)
	}
	c = cm(t, Z)
	dc = dcmdt(Z, U)
	q = c * v
}

BREAKPOINT { 
	SOLVE states METHOD derivimplicit
	i = dc * v * area * (1e-5)
}

DERIVATIVE states {
	ng' = (0.001) * sele(Z, ng)
	Z' = U
	U' = dUdt(Z, U, ng)

}

NET_RECEIVE(w) {
	LOCAL qbefore, cafter, epsilon, tt
	epsilon = 1e-10 (ms)
	if (flag == 1) { : turn on stim
		tt = tbegin
		qbefore = cm(tt - epsilon, Z)*v : charge prior to discontinuity
		cafter = cm(tt + epsilon, Z) : capacitance after discontinuity
		v = qbefore/cafter
		net_send(tdur, 2)
	} else if (flag == 2) { : turn off stim
		tt = tbegin + tdur
		qbefore = cm(tt - epsilon, Z)*v : charge prior to discontinuity
		cafter = cm(tt + epsilon, Z) : capacitance after discontinuity
	v = qbefore/cafter
	}
}
