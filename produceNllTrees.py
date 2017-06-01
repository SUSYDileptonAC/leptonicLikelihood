import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
from ROOT import gROOT, gStyle, TFile, TH1, TH1F, TH2F, TH1D, TH1I, TCanvas, TPad, TDirectory, TMath, TLorentzVector,TGraphAsymmErrors
from setTDRStyle import setTDRStyle

from messageLogger import messageLogger as log
import argparse	

import math
from math import pi, sqrt

from locations import locations
from helpers import readTrees,getFilePathsAndSampleNames

ROOT.gROOT.ProcessLine(\
					   "struct MyDileptonTreeFormat{\
						 Double_t mll;\
						 Int_t chargeProduct;\
						 Int_t nJets;\
						 Int_t nShiftedJetsJESUp;\
						 Int_t nShiftedJetsJESDown;\
						 Int_t nBJets;\
						 Int_t nBadMuonJets;\
						 Int_t nLooseLeptons;\
						 Int_t nIsoTracks;\
						 Int_t nVertices;\
						 Int_t runNr;\
						 Int_t lumiSec;\
						 ULong64_t eventNr;\
						 Double_t ht;\
						 Double_t met;\
						 Double_t metJESUp;\
						 Double_t metJESDown;\
						 Double_t caloMet;\
						 Double_t genMet;\
						 Double_t pt1;\
						 Double_t pt2;\
						 Double_t eta1;\
						 Double_t eta2;\
						 Double_t deltaPhiJetMet1;\
						 Double_t deltaPhiJetMet2;\
						 Double_t MT2;\
						 Double_t pt;\
						 Double_t sumMlb;\
						 Double_t sumMlbJESUp;\
						 Double_t sumMlbJESDown;\
						 Double_t deltaPhi;\
						 Double_t deltaR;\
						 Double_t ptPdfVal;\
						 Double_t sumMlbPdfVal;\
						 Double_t sumMlbJESUpPdfVal;\
						 Double_t sumMlbJESDownPdfVal;\
						 Double_t metPdfVal;\
						 Double_t metJESUpPdfVal;\
						 Double_t metJESDownPdfVal;\
						 Double_t deltaPhiPdfVal;\
						 Double_t nLL;\
						 Double_t nLLJESUp;\
						 Double_t nLLJESDown;\
						 Double_t weight;\
						 Double_t genWeight;\
						 Double_t bTagWeight;\
						 Double_t genPtTop1;\
						 Double_t genPtTop2;\
						 Double_t leptonFullSimScaleFactor1;\
						 Double_t leptonFullSimScaleFactor2;\
						 Double_t leptonFullSimScaleFactorErr1;\
						 Double_t leptonFullSimScaleFactorErr2;\
						 Double_t motherPdgId1;\
						 Double_t motherPdgId2;\
						 Double_t grandMotherPdgId1;\
						 Double_t grandMotherPdgId2;\
						};")

from ROOT import MyDileptonTreeFormat

def setMT2_momenta(m1,px1,py1,m2,px2,py2,metx,mety):
	global solved
	solved = False     ##reset solved tag when momenta are changed.
	global momenta_set
	momenta_set = True

	global ma
	ma = abs(m1)  ## mass cannot be negative

	if ma < 0.:
		ma = 0.

	global pax, pay, masq, Easq, Ea
	pax  = px1 
	pay  = py1
	masq = ma*ma
	Easq = masq+pax*pax+pay*pay
	Ea   = sqrt(Easq)

	global mb
	mb = abs(m2)

	if mb < 0.:
		mb = 0.

	global pbx, pby, mbsq, Ebsq, Eb
	pbx  = px2
	pby  = py2
	mbsq = mb*mb
	Ebsq = mbsq+pbx*pbx+pby*pby
	Eb   = sqrt(Ebsq)

	global pmissx, pmissy, pmissxsq,pmissysq
	pmissx = metx
	pmissy = mety
	pmissxsq = pmissx*pmissx
	pmissysq = pmissy*pmissy

## set ma>= mb
	if masq < mbsq:
		temp = pax 
		pax  = pbx 
		pbx  = temp
		
		temp = pay  
		pay  = pby 
		pby  = temp
		
		temp = Ea   
		Ea   = Eb  
		Eb   = temp
		
		temp = Easq
		Easq = Ebsq
		Ebsq = temp
		
		temp = masq
		masq = mbsq
		mbsq = temp
		
		temp = ma  
		ma   = mb 
		mb   = temp   
##normalize max{Ea, Eb} to 100
	global scale;
	if Ea > Eb:
		scale = Ea/100.
	else:
		scale = Eb/100.

	if sqrt(pmissxsq+pmissysq)/100 > scale:
		scale = sqrt(pmissxsq+pmissysq)/100

	scalesq = scale * scale
	ma  = ma/scale
	mb  = mb/scale
	masq = masq/scalesq
	mbsq = mbsq/scalesq
	pax = pax/scale
	pay = pay/scale
	pbx = pbx/scale
	pby = pby/scale
	Ea  = Ea/scale
	Eb = Eb/scale

	Easq = Easq/scalesq
	Ebsq = Ebsq/scalesq
	pmissx = pmissx/scale
	pmissy = pmissy/scale
	pmissxsq = pmissxsq/scalesq
	pmissysq = pmissysq/scalesq
	global mn
	mn   = 0.
	global mnsq
	mnsq = 0.

	global precision
	precision = 0.001


def get_mt2():
	if not momenta_set:
		print "MT2Functor::get_mt2() ==> Please set momenta first!"
		sys.exit()

	if not solved:
		mt2_bisect()
	return mt2_b*scale

def mt2_bisect():

	global solved
	solved = True

##if masses are very small, use code for massless case.  
	if masq < 0.1 and mbsq < 0.1:
		mt2_massless()
		return
	print "With mass!!!!"
	Deltasq0 = ma*(ma + 2*mn) #The minimum mass square to have two ellipses 

# find the coefficients for the two quadratic equations when Deltasq=Deltasq0.

	global mt2_b
	global a1, b1,c1, d1,e1,a2,b2,c2,d2,e2,f2
	a1 = 1-pax*pax/(Easq)
	b1 = -pax*pay/(Easq)
	c1 = 1-pay*pay/(Easq)
	d1 = -pax*(Deltasq0-masq)/(2*Easq)
	e1 = -pay*(Deltasq0-masq)/(2*Easq)
	a2 = 1-pbx*pbx/(Ebsq)
	b2 = -pbx*pby/(Ebsq)
	c2 = 1-pby*pby/(Ebsq)
	d2 = -pmissx+pbx*(Deltasq0-mbsq)/(2*Ebsq)+pbx*(pbx*pmissx+pby*pmissy)/(Ebsq)
	e2 = -pmissy+pby*(Deltasq0-mbsq)/(2*Ebsq)+pby*(pbx*pmissx+pby*pmissy)/(Ebsq)
	f2 = pmissx*pmissx+pmissy*pmissy-((Deltasq0-mbsq)/(2*Eb)+
		(pbx*pmissx+pby*pmissy)/Eb)*((Deltasq0-mbsq)/(2*Eb)+
		(pbx*pmissx+pby*pmissy)/Eb)+mnsq

# find the center of the smaller ellipse 
	x0 = (c1*d1-b1*e1)/(b1*b1-a1*c1)
	y0 = (a1*e1-b1*d1)/(b1*b1-a1*c1)


# Does the larger ellipse contain the smaller one? 
	dis=a2*x0*x0+2*b2*x0*y0+c2*y0*y0+2*d2*x0+2*e2*y0+f2

	if dis<=0.01:
		mt2_b  = sqrt(mnsq+Deltasq0)
		return


# find the coefficients for the two quadratic equations 
# coefficients for quadratic terms do not change   
# coefficients for linear and constant terms are polynomials of  
#       delta=(Deltasq-m7sq)/(2 E7sq)                               
	global d11,e11,f10,f12,d21,d20,e21,e20,f22,f21,f20
	d11 = -pax
	e11 = -pay
	f10 = mnsq
	f12 = -Easq
	d21 = (Easq*pbx)/Ebsq
	d20 = ((masq - mbsq)*pbx)/(2.*Ebsq) - pmissx + (pbx*(pbx*pmissx + pby*pmissy))/Ebsq
	e21 = (Easq*pby)/Ebsq
	e20 = ((masq - mbsq)*pby)/(2.*Ebsq) - pmissy + (pby*(pbx*pmissx + pby*pmissy))/Ebsq
	f22 = -Easq*Easq/Ebsq
	f21 = (-2*Easq*((masq - mbsq)/(2.*Eb) + (pbx*pmissx + pby*pmissy)/Eb))/Eb
	f20 = mnsq + pmissx*pmissx + pmissy*pmissy - ((masq - mbsq)/(2.*Eb) + (pbx*pmissx + pby*pmissy)/Eb)*((masq - mbsq)/(2.*Eb) + (pbx*pmissx + pby*pmissy)/Eb)

#Estimate upper bound of mT2
#when Deltasq > Deltasq_high1, the larger encloses the center of the smaller 
	p2x0 = pmissx-x0
	p2y0 = pmissy-y0
	Deltasq_high1 = 2*Eb*sqrt(p2x0*p2x0+p2y0*p2y0+mnsq)-2*pbx*p2x0-2*pby*p2y0+mbsq

#Another estimate, if both ellipses enclose the origin, Deltasq > mT2

	Deltasq_high21 = 2*Eb*sqrt(pmissx*pmissx+pmissy*pmissy+mnsq)-2*pbx*pmissx-2*pby*pmissy+mbsq
	Deltasq_high22 = 2*Ea*mn+masq

	if Deltasq_high21 < Deltasq_high22:
		Deltasq_high2 = Deltasq_high22
	else:
		Deltasq_high2 = Deltasq_high21

#pick the smaller upper bound   
	if Deltasq_high1 < Deltasq_high2:
		Deltasq_high = Deltasq_high1
	else:
		Deltasq_high = Deltasq_high2


	#lower bound
	Deltasq_low = Deltasq0;

#number of solutions at Deltasq_low should not be larger than zero
	if nsols(Deltasq_low) > 0:
# print << "MT2Functor::mt2_bisect() ==> nsolutions(Deltasq_low) > 0"
		mt2_b = sqrt(mnsq+Deltasq0)
		return

	nsols_low  = nsols(Deltasq_low)


#if nsols_high=nsols_low, we missed the region where the two ellipse overlap 
#if nsols_high=4, also need a scan because we may find the wrong tangent point.

	nsols_high = nsols(Deltasq_high)

	if nsols_high == nsols_low or nsols_high == 4:
	#foundhigh = scan_high(Deltasq_high);
		foundhigh = find_high(Deltasq_high)
		if foundhigh == 0:
			print "MT2Functor::mt2_bisect() ==> Deltasq_high not found at event "
			mt2_b = sqrt( Deltasq_low + mnsq )
			return

	while sqrt(Deltasq_high+mnsq) - sqrt(Deltasq_low+mnsq) > precision:
	#bisect
		Deltasq_mid = (Deltasq_high+Deltasq_low)/2.
		nsols_mid = nsols(Deltasq_mid)
	# if nsols_mid = 4, rescan for Deltasq_high
		if nsols_mid == 4:
			Deltasq_high = Deltasq_mid
		#scan_high(Deltasq_high);
			find_high(Deltasq_high)
			continue


		if nsols_mid != nsols_low:
			Deltasq_high = Deltasq_mid
		if nsols_mid == nsols_low:
			Deltasq_low  = Deltasq_mid
			
	mt2_b = sqrt( mnsq + Deltasq_high)
	return;

def mt2_massless():
	
	global Easq,Ebsq,Ea,Eb
	global pax,pay,pbx,pby,pmissx,pmissy
	global a2,b2,c2 
	global d21,d20,e21,e20,f22,f21,f20
	global mt2_b
	
##rotate so that pay = 0 
	theta = math.atan(pay/pax)
	s = math.sin(theta)
	c = math.cos(theta)

	
	Easq   = pax*pax+pay*pay
	Ebsq   = pbx*pbx+pby*pby
	Ea     = sqrt(Easq)
	Eb     = sqrt(Ebsq)

	pxtemp = pax*c+pay*s
	pax    = pxtemp
	pay    = 0
	pxtemp = pbx*c+pby*s
	pytemp = -s*pbx+c*pby
	pbx    = pxtemp
	pby    = pytemp
	pxtemp = pmissx*c+pmissy*s
	pytemp = -s*pmissx+c*pmissy
	pmissx = pxtemp
	pmissy = pytemp

	a2  = 1-pbx*pbx/(Ebsq)
	b2  = -pbx*pby/(Ebsq)
	c2  = 1-pby*pby/(Ebsq)

	d21 = (Easq*pbx)/Ebsq
	d20 = - pmissx +  (pbx*(pbx*pmissx + pby*pmissy))/Ebsq
	e21 = (Easq*pby)/Ebsq
	e20 = - pmissy +  (pby*(pbx*pmissx + pby*pmissy))/Ebsq
	f22 = -(Easq*Easq/Ebsq)
	f21 = -2*Easq*(pbx*pmissx + pby*pmissy)/Ebsq
	f20 = mnsq + pmissxsq + pmissysq - (pbx*pmissx + pby*pmissy)*(pbx*pmissx + pby*pmissy)/Ebsq

	Deltasq0    = 0

	Deltasq_low = Deltasq0 + precision
	nsols_low = nsols_massless(Deltasq_low)

	if nsols_low > 1: 
		mt2_b = sqrt(Deltasq0+mnsq)
		return

	#~ if( nsols_massless(Deltasq_high) > 0 )
	#~ {
		#~ mt2_b = (double) sqrt(mnsq+Deltasq0);
		#~ return;

#look for when both parablos contain origin  
	Deltasq_high1 = 2*Eb*sqrt(pmissx*pmissx+pmissy*pmissy+mnsq)-2*pbx*pmissx-2*pby*pmissy
	Deltasq_high2 = 2*Ea*mn

	if Deltasq_high1 < Deltasq_high2:
		Deltasq_high = Deltasq_high2
	else:
		Deltasq_high = Deltasq_high1

	nsols_high=nsols_massless(Deltasq_high)

	if nsols_high == nsols_low:

		foundhigh=0

		minmass  = mn
		maxmass  = sqrt(mnsq + Deltasq_high)
		mass = minmass + 0.1
		while mass < maxmass:
			Deltasq_high = mass*mass - mnsq

			nsols_high = nsols_massless(Deltasq_high)
			if nsols_high>0:
				foundhigh=1;
				Deltasq_low = (mass-0.1)*(mass-0.1) - mnsq
				break;
			mass = mass + 0.1
		if foundhigh==0:

			print "MT2Functor::mt2_massless() ==> Deltasq_high not found at event "


			mt2_b = sqrt(Deltasq_low+mnsq)
			return

	if nsols_high == nsols_low:
		print "MT2Functor::mt2_massless() ==> error: nsols_low=nsols_high= "+str(nsols_high)
		print "MT2Functor::mt2_massless() ==> Deltasq_high= "+str(Deltasq_high)
		print "MT2Functor::mt2_massless() ==> Deltasq_low= "+str(Deltasq_low << endl)

		mt2_b = sqrt(mnsq + Deltasq_low)
		return
		
	minmass = sqrt(Deltasq_low+mnsq)
	maxmass = sqrt(Deltasq_high+mnsq)
	while maxmass - minmass > precision:
		midmass   = (minmass+maxmass)/2.
		Delta_mid = midmass * midmass - mnsq
		nsols_mid = nsols_massless(Delta_mid)
		if nsols_mid != nsols_low:
			maxmass = midmass
		if nsols_mid == nsols_low:
			minmass = midmass
			
	mt2_b = minmass
	return
	
	
def nsols(Dsq):
	delta = (Dsq-masq)/(2*Easq)

#calculate coefficients for the two quadratic equations
	global d1, e1, f1, d2, e2, f2
	
	d1 = d11*delta
	e1 = e11*delta
	f1 = f12*delta*delta+f10;
	d2 = d21*delta+d20;
	e2 = e21*delta+e20;
	f2 = f22*delta*delta+f21*delta+f20;

#obtain the coefficients for the 4th order equation 
#devided by Ea^n to make the variable dimensionless

	A4 = -4*a2*b1*b2*c1 + 4*a1*b2*b2*c1 +a2*a2*c1*c1 + 4*a2*b1*b1*c2 - 4*a1*b1*b2*c2 - 2*a1*a2*c1*c2 + a1*a1*c2*c2  

	A3 =(-4*a2*b2*c1*d1 + 8*a2*b1*c2*d1 - 4*a1*b2*c2*d1 - 4*a2*b1*c1*d2 + 8*a1*b2*c1*d2 - 4*a1*b1*c2*d2 - 8*a2*b1*b2*e1 + 8*a1*b2*b2*e1 + 4*a2*a2*c1*e1 - 4*a1*a2*c2*e1 + 8*a2*b1*b1*e2 - 8*a1*b1*b2*e2 - 4*a1*a2*c1*e2 + 4*a1*a1*c2*e2)/Ea


	A2 =(4*a2*c2*d1*d1 - 4*a2*c1*d1*d2 - 4*a1*c2*d1*d2 + 4*a1*c1*d2*d2 - 8*a2*b2*d1*e1 - 8*a2*b1*d2*e1 + 16*a1*b2*d2*e1 + 4*a2*a2*e1*e1 + 16*a2*b1*d1*e2 - 8*a1*b2*d1*e2 - 8*a1*b1*d2*e2 - 8*a1*a2*e1*e2 + 4*a1*a1*e2*e2 - 4*a2*b1*b2*f1 + 4*a1*b2*b2*f1 + 2*a2*a2*c1*f1 - 2*a1*a2*c2*f1 + 4*a2*b1*b1*f2 - 4*a1*b1*b2*f2 - 2*a1*a2*c1*f2 + 2*a1*a1*c2*f2)/Easq

	A1 =(-8*a2*d1*d2*e1 + 8*a1*d2*d2*e1 + 8*a2*d1*d1*e2 - 8*a1*d1*d2*e2 - 4*a2*b2*d1*f1 - 4*a2*b1*d2*f1 + 8*a1*b2*d2*f1 + 4*a2*a2*e1*f1 - 4*a1*a2*e2*f1 + 8*a2*b1*d1*f2 - 4*a1*b2*d1*f2 - 4*a1*b1*d2*f2 -	4*a1*a2*e1*f2 + 4*a1*a1*e2*f2)/(Easq*Ea)

	A0 =(-4*a2*d1*d2*f1 + 4*a1*d2*d2*f1 + a2*a2*f1*f1 + 4*a2*d1*d1*f2 - 4*a1*d1*d2*f2 - 2*a1*a2*f1*f2 + a1*a1*f2*f2)/(Easq*Easq)

	#~ long  double A0sq, A1sq, A2sq, A3sq, A4sq;
	#~ A0sq = A0*A0;
	#~ A1sq = A1*A1;
	#~ A2sq = A2*A2;
	#~ A4sq = A4*A4;

	A3sq = A3*A3

	B3 = 4*A4
	B2 = 3*A3
	B1 = 2*A2
	B0 = A1

	C2 = -(A2/2 - 3*A3sq/(16*A4))
	C1 = -(3*A1/4. -A2*A3/(8*A4))
	C0 = -A0 + A1*A3/(16*A4)

	D1 = -B1 - (B3*C1*C1/C2 - B3*C0 -B2*C1)/C2
	D0 = -B0 - B3 *C0 *C1/(C2*C2)+ B2*C0/C2

	E0 = -C0 - C2*D0*D0/(D1*D1) + C1*D0/D1

#find the coefficients for the leading term in the Sturm sequence  
	t1 = A4
	t2 = A4
	t3 = C2
	t4 = D1
	t5 = E0


#The number of solutions depends on diffence of number of sign changes for x->Inf and x->-Inf
	nsol = signchange_n(t1,t2,t3,t4,t5) - signchange_p(t1,t2,t3,t4,t5)

#Cannot have negative number of solutions, must be roundoff effect
	if nsol < 0:
		nsol = 0

	return nsol

def nsols_massless(Dsq):
	
	#~ global d1, e1, f1, d2, e2, f2
	global d2, e2, f2

	delta = Dsq/(2*Easq)
	#~ d1    = d11*delta
	#~ e1    = e11*delta
	#~ f1    = f12*delta*delta+f10
	d2    = d21*delta+d20
	e2    = e21*delta+e20
	f2    = f22*delta*delta+f21*delta+f20
	
	if pax > 0:
		a = Ea/Dsq
	else:
		a = -Ea/Dsq
	if pax > 0:
		b = -Dsq/(4*Ea)+mnsq*Ea/Dsq
	else:
		b = Dsq/(4*Ea)-mnsq*Ea/Dsq

	A4 = a*a*a2
	A3 = 2*a*b2/Ea
	A2 = (2*a*a2*b+c2+2*a*d2)/(Easq)
	A1 = (2*b*b2+2*e2)/(Easq*Ea)
	A0 = (a2*b*b+2*b*d2+f2)/(Easq*Easq)

        A3sq = A3*A3	
        #~ long  double A0sq, A1sq, A2sq, A3sq, A4sq;
	#~ A0sq = A0*A0;
	#~ A1sq = A1*A1;
	#~ A2sq = A2*A2;
	#~ A3sq = A3*A3;
	#~ A4sq = A4*A4;

	B3 = 4*A4
	B2 = 3*A3
	B1 = 2*A2
	B0 = A1
	
	C2 = -(A2/2 - 3*A3sq/(16*A4))
	C1 = -(3*A1/4. -A2*A3/(8*A4))
	C0 = -A0 + A1*A3/(16*A4)
	
	D1 = -B1 - (B3*C1*C1/C2 - B3*C0 -B2*C1)/C2
	D0 = -B0 - B3 *C0 *C1/(C2*C2)+ B2*C0/C2

	E0 = -C0 - C2*D0*D0/(D1*D1) + C1*D0/D1

#find the coefficients for the leading term in the Sturm sequence  
	t1 = A4
	t2 = A4
	t3 = C2
	t4 = D1
	t5 = E0

	nsol = signchange_n(t1,t2,t3,t4,t5)-signchange_p(t1,t2,t3,t4,t5)
	if nsol < 0:
		nsol=0

	return nsol;


def signchange_n(t1, t2, t3, t4,t5):
	nsc=0
	if t1*t2>0:
		nsc +=1
	if t2*t3>0:
		nsc +=1
	if t3*t4>0:
		nsc +=1
	if t4*t5>0:
		nsc +=1
	return nsc;

def signchange_p(t1,t2,t3,t4,t5):
	nsc=0
	if t1*t2<0:
		nsc +=1
	if t2*t3<0:
		nsc +=1
	if t3*t4<0:
		nsc +=1
	if t4*t5<0:
		nsc +=1
	return nsc;

	
	
def find_high(Deltasq_high):
	global d1, e1, d2, e2, f2
	
	x0 = (c1*d1-b1*e1)/(b1*b1-a1*c1);
	y0 = (a1*e1-b1*d1)/(b1*b1-a1*c1);
	Deltasq_low = (mn + ma)*(mn + ma) - mnsq;
	
	while Deltasq_high - Deltasq_low > 0.001:
		Deltasq_mid = (Deltasq_high + Deltasq_low)/2.
		nsols_mid = nsols(Deltasq_mid)
		if nsols_mid == 2:
			Deltasq_high = Deltasq_mid
			return 1
		elif nsols_mid == 4:
			Deltasq_high = Deltasq_mid
			continue
		elif nsols_mid ==0:
			d1 = -pax*(Deltasq_mid-masq)/(2*Easq)
			e1 = -pay*(Deltasq_mid-masq)/(2*Easq)
			d2 = -pmissx + pbx*(Deltasq_mid - mbsq)/(2*Ebsq)+ pbx*(pbx*pmissx+pby*pmissy)/(Ebsq)
			e2 = -pmissy + pby*(Deltasq_mid - mbsq)/(2*Ebsq)+ pby*(pbx*pmissx+pby*pmissy)/(Ebsq)
			f2 = pmissx*pmissx+pmissy*pmissy-((Deltasq_mid-mbsq)/(2*Eb)+(pbx*pmissx+pby*pmissy)/Eb)*((Deltasq_mid-mbsq)/(2*Eb)+(pbx*pmissx+pby*pmissy)/Eb)+mnsq
# Does the larger ellipse contain the smaller one? 
			dis = a2*x0*x0 + 2*b2*x0*y0 + c2*y0*y0 + 2*d2*x0 + 2*e2*y0 + f2
			if dis < 0:
				Deltasq_high = Deltasq_mid
			else:
				Deltasq_low = Deltasq_mid

	return 0

	
def convertDileptonTree(tree,inputFileNames,isMC):
	pdfSourceFile = TFile(inputFileNames["pdfs"], 'READ')
	w = pdfSourceFile.Get("w")
	#~ w.Print()
	
	electronScaleFactorFile = TFile(inputFileNames["electronScaleFactorFile"])
	electronIDScaleFactorHisto = TH2F(electronScaleFactorFile.Get(inputFileNames["electronIDScaleFactorHisto"]))
	electronIsoScaleFactorHisto = TH2F(electronScaleFactorFile.Get(inputFileNames["electronIsoScaleFactorHisto"]))
	electronConvMissHitScaleFactorHisto = TH2F(electronScaleFactorFile.Get(inputFileNames["electronConvMissHitScaleFactorHisto"]))
	
	electronTrackScaleFactorFile = TFile(inputFileNames["electronTrackScaleFactorFile"])
	electronTrackScaleFactorHisto = TH2F(electronTrackScaleFactorFile.Get(inputFileNames["electronTrackScaleFactorHisto"]))
	
	muonIDScaleFactorFile = TFile(inputFileNames["muonIDScaleFactorFile"])
	muonIDScaleFactorHisto = TH2F(muonIDScaleFactorFile.Get(inputFileNames["muonIDScaleFactorHisto"]))
	
	muonIsoScaleFactorFile = TFile(inputFileNames["muonIsoScaleFactorFile"])
	muonIsoScaleFactorHisto = TH2F(muonIsoScaleFactorFile.Get(inputFileNames["muonIsoScaleFactorHisto"]))
	
	muonIP2DScaleFactorFile = TFile(inputFileNames["muonIP2DScaleFactorFile"])
	muonIP2DScaleFactorHisto = TH2F(muonIP2DScaleFactorFile.Get(inputFileNames["muonIP2DScaleFactorHisto"]))
	
	muonSIP3DScaleFactorFile = TFile(inputFileNames["muonSIP3DScaleFactorFile"])
	muonSIP3DcaleFactorHisto = TH2F(muonSIP3DScaleFactorFile.Get(inputFileNames["muonSIP3DcaleFactorHisto"]))
	
	#~ muonTrackScaleFactorFile = TFile(inputFileNames["muonTrackScaleFactorsFile"])
	#~ muonTrackScaleFactorEtaGraph = TGraphAsymmErrors(muonSIP3DScaleFactorFile.Get(inputFileNames["muonrackScaleFactorsHistoEta"]))
	#~ muonTrackScaleFactorEtaHisto = TH2F(muonSIP3DcaleFactorEtaGraph.GetHistogram())
	#~ muonTrackScaleFactorVtxGraph = TGraphAsymmErrors(muonSIP3DScaleFactorFile.Get(inputFileNames["muonrackScaleFactorsHistoNVtx"]))
	#~ muonTrackScaleFactorVtxHisto = TH2F(muonSIP3DcaleFactorVtxGraph.GetHistogram())
	
	
	### Somehow no root files of the muon track scale factors were provided and I had to read them from plotted histograms ...
	muonTrackScaleFactorEtaHisto = TH1F("muonTrackScaleFactorEtaHisto","muonTrackScaleFactorEta",12,0.,2.4)
	muonTrackScaleFactorVtxHisto = TH1F("muonTrackScaleFactorVtxHisto","muonTrackScaleFactorVtx",23,0.,46)
	
	muonTrackScaleFactorEtaHisto.SetBinContent(1,0.9970)
	muonTrackScaleFactorEtaHisto.SetBinContent(2,0.9977)
	muonTrackScaleFactorEtaHisto.SetBinContent(3,0.9981)
	muonTrackScaleFactorEtaHisto.SetBinContent(4,0.9978)
	muonTrackScaleFactorEtaHisto.SetBinContent(5,0.9980)
	muonTrackScaleFactorEtaHisto.SetBinContent(6,0.9972)
	muonTrackScaleFactorEtaHisto.SetBinContent(7,0.9962)
	muonTrackScaleFactorEtaHisto.SetBinContent(8,0.9955)
	muonTrackScaleFactorEtaHisto.SetBinContent(9,0.9958)
	muonTrackScaleFactorEtaHisto.SetBinContent(10,0.9939)
	muonTrackScaleFactorEtaHisto.SetBinContent(11,0.9929)
	muonTrackScaleFactorEtaHisto.SetBinContent(12,0.9873)
	
	muonTrackScaleFactorVtxHisto.SetBinContent(1,0.9995)
	muonTrackScaleFactorVtxHisto.SetBinContent(2,0.9993)
	muonTrackScaleFactorVtxHisto.SetBinContent(3,0.9994)
	muonTrackScaleFactorVtxHisto.SetBinContent(4,0.9991)
	muonTrackScaleFactorVtxHisto.SetBinContent(5,0.9991)
	muonTrackScaleFactorVtxHisto.SetBinContent(6,0.9985)
	muonTrackScaleFactorVtxHisto.SetBinContent(7,0.9978)
	muonTrackScaleFactorVtxHisto.SetBinContent(8,0.9972)
	muonTrackScaleFactorVtxHisto.SetBinContent(9,0.9963)
	muonTrackScaleFactorVtxHisto.SetBinContent(10,0.9959)
	muonTrackScaleFactorVtxHisto.SetBinContent(11,0.9949)
	muonTrackScaleFactorVtxHisto.SetBinContent(12,0.9941)
	muonTrackScaleFactorVtxHisto.SetBinContent(13,0.9932)
	muonTrackScaleFactorVtxHisto.SetBinContent(14,0.9929)
	muonTrackScaleFactorVtxHisto.SetBinContent(15,0.9913)
	muonTrackScaleFactorVtxHisto.SetBinContent(16,0.9902)
	muonTrackScaleFactorVtxHisto.SetBinContent(17,0.9878)
	muonTrackScaleFactorVtxHisto.SetBinContent(18,0.9877)
	muonTrackScaleFactorVtxHisto.SetBinContent(19,0.9854)
	muonTrackScaleFactorVtxHisto.SetBinContent(20,0.9882)
	muonTrackScaleFactorVtxHisto.SetBinContent(21,0.9745)
	muonTrackScaleFactorVtxHisto.SetBinContent(22,0.9676)
	muonTrackScaleFactorVtxHisto.SetBinContent(23,0.9466)
	
	# TODO: make selection more efficient
	log.logDebug("Converting DileptonTree")
	
	data = MyDileptonTreeFormat()
	newTree = ROOT.TTree("treeInvM", "Dilepton Tree")
	newTree.SetDirectory(0)
	newTree.Branch("mll", ROOT.AddressOf(data, "mll"), "mll/D")
	newTree.Branch("chargeProduct", ROOT.AddressOf(data, "chargeProduct"), "chargeProduct/I")
	newTree.Branch("nJets", ROOT.AddressOf(data, "nJets"), "nJets/I")
	newTree.Branch("nShiftedJetsJESUp", ROOT.AddressOf(data, "nShiftedJetsJESUp"), "nShiftedJetsJESUp/I")
	newTree.Branch("nShiftedJetsJESDown", ROOT.AddressOf(data, "nShiftedJetsJESDown"), "nShiftedJetsJESDown/I")
	newTree.Branch("nBJets", ROOT.AddressOf(data, "nBJets"), "nBJets/I")
	newTree.Branch("nBadMuonJets", ROOT.AddressOf(data, "nBadMuonJets"), "nBadMuonJets/I")
	newTree.Branch("nLooseLeptons", ROOT.AddressOf(data, "nLooseLeptons"), "nLooseLeptons/I")
	newTree.Branch("nIsoTracks", ROOT.AddressOf(data, "nIsoTracks"), "nIsoTracks/I")
	newTree.Branch("nVertices", ROOT.AddressOf(data, "nVertices"), "nVertices/I")
	newTree.Branch("runNr", ROOT.AddressOf(data, "runNr"), "runNr/I")
	newTree.Branch("lumiSec", ROOT.AddressOf(data, "lumiSec"), "lumiSec/I")
	newTree.Branch("eventNr", ROOT.AddressOf(data, "eventNr"), "eventNr/l")
	newTree.Branch("ht", ROOT.AddressOf(data, "ht"), "ht/D")
	newTree.Branch("met", ROOT.AddressOf(data, "met"), "met/D")
	newTree.Branch("metJESUp", ROOT.AddressOf(data, "metJESUp"), "metJESUp/D")
	newTree.Branch("metJESDown", ROOT.AddressOf(data, "metJESDown"), "metJESDown/D")
	newTree.Branch("caloMet", ROOT.AddressOf(data, "caloMet"), "caloMet/D")
	newTree.Branch("genMet", ROOT.AddressOf(data, "genMet"), "genMet/D")
	newTree.Branch("pt1", ROOT.AddressOf(data, "pt1"), "pt1/D")
	newTree.Branch("pt2", ROOT.AddressOf(data, "pt2"), "pt2/D")
	newTree.Branch("eta1", ROOT.AddressOf(data, "eta1"), "eta1/D")
	newTree.Branch("eta2", ROOT.AddressOf(data, "eta2"), "eta2/D")
	newTree.Branch("deltaPhiJetMet1", ROOT.AddressOf(data, "deltaPhiJetMet1"), "deltaPhiJetMet1/D")
	newTree.Branch("deltaPhiJetMet2", ROOT.AddressOf(data, "deltaPhiJetMet2"), "deltaPhiJetMet2/D")
	newTree.Branch("MT2", ROOT.AddressOf(data, "MT2"), "MT2/D")
	newTree.Branch("pt", ROOT.AddressOf(data, "pt"), "pt/D")
	newTree.Branch("sumMlb", ROOT.AddressOf(data, "sumMlb"), "sumMlb/D")
	newTree.Branch("sumMlbJESUp", ROOT.AddressOf(data, "sumMlbJESUp"), "sumMlbJESUp/D")
	newTree.Branch("sumMlbJESDown", ROOT.AddressOf(data, "sumMlbJESDown"), "sumMlbJESDown/D")
	newTree.Branch("deltaPhi", ROOT.AddressOf(data, "deltaPhi"), "deltaPhi/D")
	newTree.Branch("deltaR", ROOT.AddressOf(data, "deltaR"), "deltaR/D")
	newTree.Branch("ptPdfVal", ROOT.AddressOf(data, "ptPdfVal"), "ptPdfVal/D")
	newTree.Branch("sumMlbPdfVal", ROOT.AddressOf(data, "sumMlbPdfVal"), "sumMlbPdfVal/D")
	newTree.Branch("sumMlbJESUpPdfVal", ROOT.AddressOf(data, "sumMlbJESUpPdfVal"), "sumMlbJESUpPdfVal/D")
	newTree.Branch("sumMlbJESDownPdfVal", ROOT.AddressOf(data, "sumMlbJESDownPdfVal"), "sumMlbJESDownPdfVal/D")
	newTree.Branch("deltaPhiPdfVal", ROOT.AddressOf(data, "deltaPhiPdfVal"), "deltaPhiPdfVal/D")
	newTree.Branch("metPdfVal", ROOT.AddressOf(data, "metPdfVal"), "metPdfVal/D")
	newTree.Branch("metJESUpPdfVal", ROOT.AddressOf(data, "metJESUpPdfVal"), "metJESUpPdfVal/D")
	newTree.Branch("metJESDownPdfVal", ROOT.AddressOf(data, "metJESDownPdfVal"), "metJESDownPdfVal/D")
	newTree.Branch("nLL", ROOT.AddressOf(data, "nLL"), "nLL/D")
	newTree.Branch("nLLJESUp", ROOT.AddressOf(data, "nLLJESUp"), "nLLJESUp/D")
	newTree.Branch("nLLJESDown", ROOT.AddressOf(data, "nLLJESDown"), "nLLJESDown/D")
	newTree.Branch("weight", ROOT.AddressOf(data, "weight"), "weight/D")
	#~ newTree.Branch("weightUp", ROOT.AddressOf(data, "weightUp"), "weightUp/D")
	#~ newTree.Branch("weightDown", ROOT.AddressOf(data, "weightDown"), "weightDown/D")
	newTree.Branch("genWeight", ROOT.AddressOf(data, "genWeight"), "genWeight/D")
	newTree.Branch("bTagWeight", ROOT.AddressOf(data, "bTagWeight"), "bTagWeight/D")
	newTree.Branch("genPtTop1", ROOT.AddressOf(data, "genPtTop1"), "genPtTop1/D")
	newTree.Branch("genPtTop2", ROOT.AddressOf(data, "genPtTop2"), "genPtTop2/D")
	newTree.Branch("leptonFullSimScaleFactor1", ROOT.AddressOf(data, "leptonFullSimScaleFactor1"), "leptonFullSimScaleFactor1/D")
	newTree.Branch("leptonFullSimScaleFactor2", ROOT.AddressOf(data, "leptonFullSimScaleFactor2"), "leptonFullSimScaleFactor2/D")
	newTree.Branch("leptonFullSimScaleFactorErr1", ROOT.AddressOf(data, "leptonFullSimScaleFactorErr1"), "leptonFullSimScaleFactorErr1/D")
	newTree.Branch("leptonFullSimScaleFactorErr2", ROOT.AddressOf(data, "leptonFullSimScaleFactorErr2"), "leptonFullSimScaleFactorErr2/D")
	newTree.Branch("motherPdgId1", ROOT.AddressOf(data, "motherPdgId1"), "motherPdgId1/D")
	newTree.Branch("motherPdgId2", ROOT.AddressOf(data, "motherPdgId2"), "motherPdgId2/D")
	newTree.Branch("grandMotherPdgId1", ROOT.AddressOf(data, "grandMotherPdgId1"), "grandMotherPdgId1/D")
	newTree.Branch("grandMotherPdgId2", ROOT.AddressOf(data, "grandMotherPdgId2"), "grandMotherPdgId2/D")


	# only part of tree?
	iMax = tree.GetEntries()
	
	# Fill tree
	for i in xrange(iMax):
		if (tree.GetEntry(i) > 0):
			data.weight = tree.weight
			data.genWeight = tree.genWeight
			if isMC:
				data.bTagWeight = tree.bTagWeight
			else:
				data.bTagWeight = 1.
			data.genPtTop1 = tree.genPtTop1
			data.genPtTop2 = tree.genPtTop2
			if not isMC:
				data.leptonFullSimScaleFactor1 = 1.
				data.leptonFullSimScaleFactor2 = 1.
				data.leptonFullSimScaleFactorErr1 = 0.
				data.leptonFullSimScaleFactorErr2 = 0.
			else:
				fullSimSF1 = 1.
				fullSimSF2 = 1.
				fullSimSFErr1 = 0.
				fullSimSFErr2 = 0.
				
				if abs(tree.pdgId1) == 11:
					if tree.pt1 > 200.:
						tempPt = 199.
					else:
						tempPt = tree.pt1
					if tree.pt1 > 500.:
						tempPtTrack = 499.
					elif tree.pt1 < 25.:
						tempPtTrack = 26.
					else:
						tempPtTrack = tree.pt1
					tempEta = abs(tree.eta1)
					
					fullSimSF1 = fullSimSF1 * electronIDScaleFactorHisto.GetBinContent(electronIDScaleFactorHisto.GetXaxis().FindBin(tempPt),electronIDScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF1 = fullSimSF1 * electronIsoScaleFactorHisto.GetBinContent(electronIsoScaleFactorHisto.GetXaxis().FindBin(tempPt),electronIsoScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF1 = fullSimSF1 * electronConvMissHitScaleFactorHisto.GetBinContent(electronConvMissHitScaleFactorHisto.GetXaxis().FindBin(tempPt),electronConvMissHitScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF1 = fullSimSF1 * electronTrackScaleFactorHisto.GetBinContent(electronTrackScaleFactorHisto.GetXaxis().FindBin(tree.eta1),electronTrackScaleFactorHisto.GetYaxis().FindBin(tempPtTrack))
					
					
					## get the uncertainties from the SF histograms for ee + 1% for the track reconstruction SF
					fullSimSFErr1 = electronIDScaleFactorHisto.GetBinError(electronIDScaleFactorHisto.GetXaxis().FindBin(tempPt),electronIDScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSFErr1 = sqrt(fullSimSFErr1**2 + electronIsoScaleFactorHisto.GetBinError(electronIsoScaleFactorHisto.GetXaxis().FindBin(tempPt),electronIsoScaleFactorHisto.GetYaxis().FindBin(tempEta))**2)
					fullSimSFErr1 = sqrt(fullSimSFErr1**2 + electronConvMissHitScaleFactorHisto.GetBinError(electronConvMissHitScaleFactorHisto.GetXaxis().FindBin(tempPt),electronConvMissHitScaleFactorHisto.GetYaxis().FindBin(tempEta))**2)
					
					fullSimSFErr1 =  sqrt(fullSimSFErr1**2 + 0.01**2)
					
				if abs(tree.pdgId1) == 13:
					if tree.pt1 > 120.:
						tempPt = 119.
					else:
						tempPt = tree.pt1
					#~ if tree.pt1 > 500.:
						#~ tempPtTrack = 199.
					#~ else:
						#~ tempPtTrack = tree.pt1
					tempEta = abs(tree.eta1)
					if tree.nVertices > 45:
						tempnVtx = 45
					else:
						tempnVtx = tree.nVertices
					
					fullSimSF1 = fullSimSF1 * muonIDScaleFactorHisto.GetBinContent(muonIDScaleFactorHisto.GetXaxis().FindBin(tempPt),muonIDScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF1 = fullSimSF1 * muonIsoScaleFactorHisto.GetBinContent(muonIsoScaleFactorHisto.GetXaxis().FindBin(tempPt),muonIsoScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF1 = fullSimSF1 * muonIP2DScaleFactorHisto.GetBinContent(muonIP2DScaleFactorHisto.GetXaxis().FindBin(tempPt),muonIP2DScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF1 = fullSimSF1 * muonSIP3DcaleFactorHisto.GetBinContent(muonSIP3DcaleFactorHisto.GetXaxis().FindBin(tempPt),muonSIP3DcaleFactorHisto.GetYaxis().FindBin(tempEta))
					
					fullSimSF1 = fullSimSF1 * muonTrackScaleFactorEtaHisto.GetBinContent(muonTrackScaleFactorEtaHisto.GetYaxis().FindBin(tempEta))
					fullSimSF1 = fullSimSF1 * muonTrackScaleFactorVtxHisto.GetBinContent(muonTrackScaleFactorVtxHisto.GetYaxis().FindBin(tempnVtx))
					
					### constant 3% uncertainty for muons
					fullSimSFErr1 = 0.03
				
				if abs(tree.pdgId2) == 11:
					if tree.pt2 > 200.:
						tempPt = 199.
					else:
						tempPt = tree.pt2
					if tree.pt2 > 500.:
						tempPtTrack = 499.
					elif tree.pt2 < 25.:
						tempPtTrack = 26.
					else:
						tempPtTrack = tree.pt2
					tempEta = abs(tree.eta2)
					
					fullSimSF2 = fullSimSF2 * electronIDScaleFactorHisto.GetBinContent(electronIDScaleFactorHisto.GetXaxis().FindBin(tempPt),electronIDScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF2 = fullSimSF2 * electronIsoScaleFactorHisto.GetBinContent(electronIsoScaleFactorHisto.GetXaxis().FindBin(tempPt),electronIsoScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF2 = fullSimSF2 * electronConvMissHitScaleFactorHisto.GetBinContent(electronConvMissHitScaleFactorHisto.GetXaxis().FindBin(tempPt),electronConvMissHitScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF2 = fullSimSF2 * electronTrackScaleFactorHisto.GetBinContent(electronTrackScaleFactorHisto.GetXaxis().FindBin(tree.eta2),electronTrackScaleFactorHisto.GetYaxis().FindBin(tempPtTrack))
					
					## get the uncertainties from the SF histograms for ee + 1% for the track reconstruction SF
					fullSimSFErr2 = electronIDScaleFactorHisto.GetBinError(electronIDScaleFactorHisto.GetXaxis().FindBin(tempPt),electronIDScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSFErr2 = sqrt(fullSimSFErr2**2 + electronIsoScaleFactorHisto.GetBinError(electronIsoScaleFactorHisto.GetXaxis().FindBin(tempPt),electronIsoScaleFactorHisto.GetYaxis().FindBin(tempEta))**2)
					fullSimSFErr2 = sqrt(fullSimSFErr2**2 + electronConvMissHitScaleFactorHisto.GetBinError(electronConvMissHitScaleFactorHisto.GetXaxis().FindBin(tempPt),electronConvMissHitScaleFactorHisto.GetYaxis().FindBin(tempEta))**2)
					
					fullSimSFErr2 =  sqrt(fullSimSFErr2**2 + 0.01**2)
					
				if abs(tree.pdgId2) == 13:
					if tree.pt2 > 120.:
						tempPt = 119.
					else:
						tempPt = tree.pt2
					#~ if tree.pt2 > 500.:
						#~ tempPtTrack = 199.
					#~ else:
						#~ tempPtTrack = tree.pt2
					tempEta = abs(tree.eta2)
					
					if tree.nVertices > 45:
						tempnVtx = 45
					else:
						tempnVtx = tree.nVertices
					
					fullSimSF2 = fullSimSF2 * muonIDScaleFactorHisto.GetBinContent(muonIDScaleFactorHisto.GetXaxis().FindBin(tempPt),muonIDScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF2 = fullSimSF2 * muonIsoScaleFactorHisto.GetBinContent(muonIsoScaleFactorHisto.GetXaxis().FindBin(tempPt),muonIsoScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF2 = fullSimSF2 * muonIP2DScaleFactorHisto.GetBinContent(muonIP2DScaleFactorHisto.GetXaxis().FindBin(tempPt),muonIP2DScaleFactorHisto.GetYaxis().FindBin(tempEta))
					fullSimSF2 = fullSimSF2 * muonSIP3DcaleFactorHisto.GetBinContent(muonSIP3DcaleFactorHisto.GetXaxis().FindBin(tempPt),muonSIP3DcaleFactorHisto.GetYaxis().FindBin(tempEta))
					
					fullSimSF2 = fullSimSF2 * muonTrackScaleFactorEtaHisto.GetBinContent(muonTrackScaleFactorEtaHisto.GetYaxis().FindBin(tempEta))
					fullSimSF2 = fullSimSF2 * muonTrackScaleFactorVtxHisto.GetBinContent(muonTrackScaleFactorVtxHisto.GetYaxis().FindBin(tempnVtx))
					
					### constant 3% uncertainty for muons
					fullSimSFErr2 = 0.03
					
				data.leptonFullSimScaleFactor1 = fullSimSF1
				data.leptonFullSimScaleFactor2 = fullSimSF2
				data.leptonFullSimScaleFactorErr1 = fullSimSFErr1
				data.leptonFullSimScaleFactorErr2 = fullSimSFErr2
			data.motherPdgId1 = tree.motherPdgId1
			data.motherPdgId2 = tree.motherPdgId2
			data.mll = tree.p4.M()
			data.chargeProduct = tree.chargeProduct
			data.nJets = tree.nJets
			data.nShiftedJetsJESUp = tree.nShiftedJetsJESUp
			data.nShiftedJetsJESDown = tree.nShiftedJetsJESDown
			data.nBJets = tree.nBJets
			data.nBadMuonJets = tree.nBadMuonJets
			data.nLooseLeptons = tree.nLooseLeptons
			data.nIsoTracks = tree.nIsoTracks
			data.nVertices = tree.nVertices
			data.runNr = tree.runNr
			data.lumiSec = tree.lumiSec
			data.eventNr = tree.eventNr
			data.ht = tree.ht
			data.met = tree.met
			data.metJESUp = tree.metJESUp
			data.metJESDown = tree.metJESDown
			data.caloMet = tree.caloMet
			data.genMet = tree.genMet
			data.pt1 = tree.pt1
			data.pt2 = tree.pt2
			data.eta1 = tree.eta1
			data.eta2 = tree.eta2
			data.deltaPhiJetMet1 = tree.deltaPhiJetMet1
			data.deltaPhiJetMet2 = tree.deltaPhiJetMet2
			### In some cases the MT2 calculation in the tree maker screwed things up and produced inf or nan. Recalculate in such cases 
			### (should be fixed now, but does not hurt to ask if MT2 was calculated correctly)
			if math.isnan(tree.MT2) or math.isinf(tree.MT2) or tree.MT2 > 10000.:
				setMT2_momenta(tree.lepton1.M(),tree.lepton1.Px(),tree.lepton1.Py(),tree.lepton2.M(),tree.lepton2.Px(),tree.lepton2.Py(),tree.vMet.Px(),tree.vMet.Py())
				MT2Val = get_mt2()
				data.MT2 = MT2Val
			else:
				data.MT2 = tree.MT2
			data.pt = tree.p4.Pt()
			data.sumMlb = tree.sumMlb
			data.sumMlbJESUp = tree.sumMlbJESUp
			data.sumMlbJESDown = tree.sumMlbJESDown
			data.deltaPhi = abs(tree.deltaPhi)
			data.deltaR = tree.deltaR
			if isMC:					
				metPdf = w.pdf("met_analyticalPDF_MC")
				ptPdf = w.pdf("zpt_analyticalPDF_MC")
				sumMlbPdf = w.pdf("mlb_analyticalPDF_MC")
				deltaPhiPdf = w.pdf("ldp_analyticalPDF_MC")
			else:
				metPdf = w.pdf("met_analyticalPDF_DA")
				ptPdf = w.pdf("zpt_analyticalPDF_DA")
				sumMlbPdf = w.pdf("mlb_analyticalPDF_DA")
				deltaPhiPdf = w.pdf("ldp_analyticalPDF_DA")
			w.var("lepsZPt_Edge").setVal(tree.p4.Pt())
			w.var("sum_mlb_Edge").setVal(tree.sumMlb)
			w.var("lepsDPhi_Edge").setVal(abs(tree.deltaPhi))
			w.var("met_Edge").setVal(tree.met)
			obsPt = ROOT.RooArgSet(w.var("lepsZPt_Edge"))
			obsMlb = ROOT.RooArgSet(w.var("sum_mlb_Edge"))
			obsDeltaPhi = ROOT.RooArgSet(w.var("lepsDPhi_Edge"))
			obsMet = ROOT.RooArgSet(w.var("met_Edge"))
			data.ptPdfVal = ptPdf.getVal(obsPt)
			data.sumMlbPdfVal = sumMlbPdf.getVal(obsMlb)
			data.deltaPhiPdfVal = deltaPhiPdf.getVal(obsDeltaPhi)
			data.metPdfVal = metPdf.getVal(obsMet)
			
			w.var("met_Edge").setVal(tree.metJESUp)
			w.var("sum_mlb_Edge").setVal(tree.sumMlbJESUp)	
			obsMet = ROOT.RooArgSet(w.var("met_Edge"))
			obsMlb = ROOT.RooArgSet(w.var("sum_mlb_Edge"))
			data.metJESUpPdfVal = metPdf.getVal(obsMet)
			data.sumMlbJESUpPdfVal = sumMlbPdf.getVal(obsMlb)	
			w.var("met_Edge").setVal(tree.metJESDown)
			w.var("sum_mlb_Edge").setVal(tree.sumMlbJESDown)
			obsMet = ROOT.RooArgSet(w.var("met_Edge"))
			obsMlb = ROOT.RooArgSet(w.var("sum_mlb_Edge"))
			data.metJESDownPdfVal = metPdf.getVal(obsMet)	
			data.sumMlbJESDownPdfVal = sumMlbPdf.getVal(obsMlb)	
			
			
			data.nLL = - TMath.Log(data.ptPdfVal) - TMath.Log(data.sumMlbPdfVal) - TMath.Log(data.metPdfVal) - TMath.Log(data.deltaPhiPdfVal)
			data.nLLJESUp = - TMath.Log(data.ptPdfVal) - TMath.Log(data.sumMlbJESUpPdfVal) - TMath.Log(data.metJESUpPdfVal) - TMath.Log(data.deltaPhiPdfVal)
			data.nLLJESDown = - TMath.Log(data.ptPdfVal) - TMath.Log(data.sumMlbJESDownPdfVal) - TMath.Log(data.metJESDownPdfVal) - TMath.Log(data.deltaPhiPdfVal)
			newTree.Fill()

	return newTree

					

def main():
	
	
	

	parser = argparse.ArgumentParser(description='edge fitter reloaded.')
	
	parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False,
						  help="Verbose mode.")
	parser.add_argument("-S", "--sample", action="store", dest="sample", default="all",
						  help="choose a sample to convert, default is all.")
	parser.add_argument("-s", "--selection", dest = "selection" , action="store", default="SignalInclusive",
						  help="selection which to apply.")	
					
	args = parser.parse_args()	
	
	canv = TCanvas("canv", "canv",800,800)
	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
	style=setTDRStyle()	
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()
	
	cut = " chargeProduct < 0 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25))  && abs(eta1)<2.4  && abs(eta2) < 2.4 && ((abs(eta1) < 1.4 || abs(eta1) > 1.6) && (abs(eta2) < 1.4 || abs(eta2) > 1.6)) && deltaR > 0.1 && p4.M() > 20"
	
	### Likelihood does actually only make sense for met > 150, use met > 100 to be able to do the RSFOF calculation on the resulting 
	### smaller samples as well which is much faster. Remember not to use the likelihood for met < 150
	cut += " && nJets >= 2 && met > 100"
	#~ cut = ""
	
	



	region = args.selection
	
	from defs import Regions
	if not region in dir(Regions):
		print "invalid region, exiting"
		sys.exit()
	else:	
		Selection = getattr(Regions,region)
		
	path = locations.dataSetPath
	outpath = locations.dataSetPathNLL
	swversion = "sw8026v1001"
	cutsversion = "cutsV34"
	inputFiles = {
					"pdfs":"pdfs_forMoriond_ver5.root",
					
					"electronScaleFactorFile":"ScaleFactorsElectrons.root",
					"electronIDScaleFactorHisto":"GsfElectronToMVATightTightIP2DSIP3D4",
					"electronIsoScaleFactorHisto":"MVAVLooseElectronToMini",
					"electronConvMissHitScaleFactorHisto":"MVATightElectronToConvVetoIHit0",
					
					"electronTrackScaleFactorFile":"TrackScaleFactorsElectrons.root",
					"electronTrackScaleFactorHisto":"EGamma_SF2D",
	
					"muonIDScaleFactorFile":"ScaleFactorMuonID.root",
					"muonIDScaleFactorHisto":"SF",
					"muonIsoScaleFactorFile":"ScaleFactorMuonMiniIso.root",
					"muonIsoScaleFactorHisto":"SF",
					"muonIP2DScaleFactorFile":"ScaleFactorMuonIP2D.root",
					"muonIP2DScaleFactorHisto":"SF",
					"muonSIP3DScaleFactorFile":"ScaleFactorMuonSIP3D.root",
					"muonSIP3DcaleFactorHisto":"SF",
					
					"muonTrackScaleFactorsFile":"TrackScaleFactorsMuons.root",
					"muonrackScaleFactorsHistoEta":"ratio_eff_eta3_dr030e030_corr",
					"muonrackScaleFactorsHistoNVtx":"ratio_eff_vtx_dr030e030_corr",
					}
	
	
	# init ROOT
	gROOT.Reset()
	
	samples = getFilePathsAndSampleNames(path)
	
	treeEEraw = readTrees(path, "EE")
	treeEMuraw = readTrees(path, "EMu")
	treeMuMuraw = readTrees(path, "MuMu")
	
	### Only use events where the MET filters are ok
	baseCut = cut+ " && Flag_HBHENoiseFilter > 0 && Flag_HBHENoiseIsoFilter > 0 && Flag_globalTightHalo2016Filter > 0 && Flag_goodVertices > 0 && Flag_EcalDeadCellTriggerPrimitiveFilter > 0"
	
	for sample in samples:
		if sample == "WWTo4Q_Powheg_Spring16_25ns" or sample == "WWTo4Q_Powheg_Summer16_25ns":  ### No events remain
			continue
		if sample == "MergedData":
			#~ continue
			isMC = False
			cut = baseCut+ " && Flag_eeBadScFilter > 0" ### additional MET filter for data
		else:
			isMC = True
			cut = baseCut
			#~ if not sample in ["WZTo3LNu_Powheg_Summer16_25ns","ZZTo2L2Nu_Powheg_Summer16_25ns","TTZToLLNuNu_aMCatNLO_FXFX_Summer16_25ns","ST_top_tWllChannel_5f_Powheg_MadGraph_Summer16_25ns","WWZ_aMCatNLO_FXFX_Summer16_25ns","WZZ_aMCatNLO_FXFX_Summer16_25ns","ZZZ_aMCatNLO_FXFX_Summer16_25ns","TZQ_LL_aMCatNLO_Summer16_25ns","ZZTo4L_Powheg_Summer16_25ns","ZZTo2L2Q_aMCatNLO_Summer16_25ns","ZJets_Madgraph_Summer16_25ns","AStar_Madgraph_Summer16_25ns"]:
			#~ if not sample in ["ST_top_tWllChannel_5f_Powheg_MadGraph_Summer16_25ns","TZQ_LL_aMCatNLO_Summer16_25ns","WZZ_aMCatNLO_FXFX_Summer16_25ns","WWZ_aMCatNLO_FXFX_Summer16_25ns","ZZZ_aMCatNLO_FXFX_Summer16_25ns","TTZToLLNuNu_aMCatNLO_FXFX_Summer16_25ns","ZZTo2L2Nu_Powheg_Summer16_25ns","WZTo3LNu_Powheg_Summer16_25ns"]:
				#~ continue
			
			
		for sampleName, filePath in getFilePathsAndSampleNames(path).iteritems():
			if sampleName == sample:
				inFile = TFile(filePath, "read")
				denominator = inFile.FindObjectAny("analysis paths").GetBinContent(1)
				histo = TH1F("analysis paths", "analysis paths", 2, 0 ,1)
				histo.GetXaxis().SetBinLabel(1, "None")
				histo.SetBinContent(1, denominator)
				 
		
		#~ if sample == "MergedData":
		print sample
		f1 = TFile("%s/%s.processed.%s.root"%(outpath,swversion,sample),"RECREATE")
		folder = f1.mkdir("%sDileptonFinalTrees"%cutsversion)
		
		### Only use events with trigger fired
		cutEE = cut+" && (HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v > 0 || HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v > 0)"
		cutEM = cut+" && (HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v > 0 || HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v > 0)"
		cutMM = cut+" && (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v > 0 || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v > 0 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v > 0 || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v > 0 || HLT_Mu27_TkMu8_v > 0 || HLT_Mu30_TkMu11_v > 0)"

			
			
		treeEECut = treeEEraw[sample].CopyTree(cutEE)
		treeEMuCut = treeEMuraw[sample].CopyTree(cutEM)
		treeMuMuCut = treeMuMuraw[sample].CopyTree(cutMM)
		
		
		treeEEConverted = convertDileptonTree(treeEECut,inputFiles,isMC)
		treeEMuConverted = convertDileptonTree(treeEMuCut,inputFiles,isMC)
		treeMuMuConverted = convertDileptonTree(treeMuMuCut,inputFiles,isMC)
		
		treeEEConverted.SetName("EEDileptonTree")
		treeEMuConverted.SetName("EMuDileptonTree")
		treeMuMuConverted.SetName("MuMuDileptonTree")

		
		
		folder.cd()
		
		treeEEConverted.Write()
		treeEMuConverted.Write()
		treeMuMuConverted.Write()
		
		folder2 = f1.mkdir("%sDileptonCounters"%cutsversion)
		folder2.cd()
		histo.Write()
		
				
				
				
				
		f1.Close()		
		
			
			
	
			

main()
