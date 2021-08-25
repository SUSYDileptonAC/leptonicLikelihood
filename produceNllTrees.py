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
from helpers import readTrees,getFilePathsAndSampleNames,ensurePathExists

ROOT.gROOT.ProcessLine(\
                                           "struct MyDileptonTreeFormat{\
                                                 Float_t mll;\
                                                 Float_t chargeProduct;\
                                                 Int_t vetoHEM;\
                                                 Int_t affectedHEM;\
                                                 Int_t nJets;\
                                                 Int_t nShiftedJetsJESUp;\
                                                 Int_t nShiftedJetsJESDown;\
                                                 Int_t nBJets;\
                                                 Int_t nBadMuonJets;\
                                                 Int_t nLooseLeptons;\
                                                 Int_t nVertices;\
                                                 Int_t runNr;\
                                                 Int_t lumiSec;\
                                                 ULong64_t eventNr;\
                                                 Float_t ht;\
                                                 Float_t met;\
                                                 Float_t metJESUp;\
                                                 Float_t metJESDown;\
                                                 Float_t caloMet;\
                                                 Float_t genMet;\
                                                 Float_t pt1;\
                                                 Float_t pt2;\
                                                 Float_t eta1;\
                                                 Float_t eta2;\
                                                 Float_t deltaPhiJetMet1;\
                                                 Float_t deltaPhiJetMet2;\
                                                 Float_t MT2;\
                                                 Float_t pt;\
                                                 Float_t sumMlb;\
                                                 Float_t sumMlbJESUp;\
                                                 Float_t sumMlbJESDown;\
                                                 Float_t deltaPhi;\
                                                 Float_t deltaR;\
                                                 Float_t ptPdfVal;\
                                                 Float_t sumMlbPdfVal;\
                                                 Float_t sumMlbJESUpPdfVal;\
                                                 Float_t sumMlbJESDownPdfVal;\
                                                 Float_t metPdfVal;\
                                                 Float_t metJESUpPdfVal;\
                                                 Float_t metJESDownPdfVal;\
                                                 Float_t deltaPhiPdfVal;\
                                                 Float_t nLL;\
                                                 Float_t nLLJESUp;\
                                                 Float_t nLLJESDown;\
                                                 Float_t weight;\
                                                 Float_t prefireWeight;\
                                                 Float_t genWeight;\
                                                 Float_t bTagWeight;\
                                                 Float_t bTagWeightLoose;\
                                                 Float_t genPtTop1;\
                                                 Float_t genPtTop2;\
                                                 Float_t leptonFullSimScaleFactor1;\
                                                 Float_t leptonFullSimScaleFactor2;\
                                                 Float_t leptonFullSimScaleFactorErr1;\
                                                 Float_t leptonFullSimScaleFactorErr2;\
                                                 Float_t kfactor_zzto4l_pt;\
                                                 Float_t kfactor_zzto2l_pt;\
                                                 Float_t kfactor_zzto4l_mass;\
                                                 Float_t kfactor_zzto2l_mass;\
                                                 Int_t motherPdgId1;\
                                                 Int_t motherPdgId2;\
                                                 Int_t grandMotherPdgId1;\
                                                 Int_t grandMotherPdgId2;\
                                                };") # Int_t nIsoTracks;\

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

        A1 =(-8*a2*d1*d2*e1 + 8*a1*d2*d2*e1 + 8*a2*d1*d1*e2 - 8*a1*d1*d2*e2 - 4*a2*b2*d1*f1 - 4*a2*b1*d2*f1 + 8*a1*b2*d2*f1 + 4*a2*a2*e1*f1 - 4*a1*a2*e2*f1 + 8*a2*b1*d1*f2 - 4*a1*b2*d1*f2 - 4*a1*b1*d2*f2 -   4*a1*a2*e1*f2 + 4*a1*a1*e2*f2)/(Easq*Ea)

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

data = MyDileptonTreeFormat()
def convertDileptonTree(tree,inputFileNames,isMC,year):
        
        pdfSourceFile = TFile(inputFileNames["pdfs"][year], 'READ')
        w = pdfSourceFile.Get("w")
                
        # TODO: make selection more efficient
        log.logDebug("Converting DileptonTree")
        
        global data
        newTree = ROOT.TTree("treeInvM", "Dilepton Tree")
        newTree.SetDirectory(0)
        newTree.Branch("mll", ROOT.AddressOf(data, "mll"), "mll/F")
        newTree.Branch("chargeProduct", ROOT.AddressOf(data, "chargeProduct"), "chargeProduct/F")
        newTree.Branch("vetoHEM", ROOT.AddressOf(data, "vetoHEM"), "vetoHEM/I")
        newTree.Branch("affectedHEM", ROOT.AddressOf(data, "affectedHEM"), "affectedHEM/I")
        newTree.Branch("nJets", ROOT.AddressOf(data, "nJets"), "nJets/I")
        newTree.Branch("nShiftedJetsJESUp", ROOT.AddressOf(data, "nShiftedJetsJESUp"), "nShiftedJetsJESUp/I")
        newTree.Branch("nShiftedJetsJESDown", ROOT.AddressOf(data, "nShiftedJetsJESDown"), "nShiftedJetsJESDown/I")
        newTree.Branch("nBJets", ROOT.AddressOf(data, "nBJets"), "nBJets/I")
        newTree.Branch("nBadMuonJets", ROOT.AddressOf(data, "nBadMuonJets"), "nBadMuonJets/I")
        newTree.Branch("nLooseLeptons", ROOT.AddressOf(data, "nLooseLeptons"), "nLooseLeptons/I")
        #newTree.Branch("nIsoTracks", ROOT.AddressOf(data, "nIsoTracks"), "nIsoTracks/I")
        newTree.Branch("nVertices", ROOT.AddressOf(data, "nVertices"), "nVertices/I")
        newTree.Branch("runNr", ROOT.AddressOf(data, "runNr"), "runNr/I")
        newTree.Branch("lumiSec", ROOT.AddressOf(data, "lumiSec"), "lumiSec/I")
        newTree.Branch("eventNr", ROOT.AddressOf(data, "eventNr"), "eventNr/l")
        newTree.Branch("ht", ROOT.AddressOf(data, "ht"), "ht/F")
        newTree.Branch("met", ROOT.AddressOf(data, "met"), "met/F")
        newTree.Branch("metJESUp", ROOT.AddressOf(data, "metJESUp"), "metJESUp/F")
        newTree.Branch("metJESDown", ROOT.AddressOf(data, "metJESDown"), "metJESDown/F")
        newTree.Branch("caloMet", ROOT.AddressOf(data, "caloMet"), "caloMet/F")
        newTree.Branch("genMet", ROOT.AddressOf(data, "genMet"), "genMet/F")
        newTree.Branch("pt1", ROOT.AddressOf(data, "pt1"), "pt1/F")
        newTree.Branch("pt2", ROOT.AddressOf(data, "pt2"), "pt2/F")
        newTree.Branch("eta1", ROOT.AddressOf(data, "eta1"), "eta1/F")
        newTree.Branch("eta2", ROOT.AddressOf(data, "eta2"), "eta2/F")
        newTree.Branch("deltaPhiJetMet1", ROOT.AddressOf(data, "deltaPhiJetMet1"), "deltaPhiJetMet1/F")
        newTree.Branch("deltaPhiJetMet2", ROOT.AddressOf(data, "deltaPhiJetMet2"), "deltaPhiJetMet2/F")
        newTree.Branch("MT2", ROOT.AddressOf(data, "MT2"), "MT2/F")
        newTree.Branch("pt", ROOT.AddressOf(data, "pt"), "pt/F")
        newTree.Branch("sumMlb", ROOT.AddressOf(data, "sumMlb"), "sumMlb/F")
        newTree.Branch("sumMlbJESUp", ROOT.AddressOf(data, "sumMlbJESUp"), "sumMlbJESUp/F")
        newTree.Branch("sumMlbJESDown", ROOT.AddressOf(data, "sumMlbJESDown"), "sumMlbJESDown/F")
        newTree.Branch("deltaPhi", ROOT.AddressOf(data, "deltaPhi"), "deltaPhi/F")
        newTree.Branch("deltaR", ROOT.AddressOf(data, "deltaR"), "deltaR/F")
        newTree.Branch("ptPdfVal", ROOT.AddressOf(data, "ptPdfVal"), "ptPdfVal/F")
        newTree.Branch("sumMlbPdfVal", ROOT.AddressOf(data, "sumMlbPdfVal"), "sumMlbPdfVal/F")
        newTree.Branch("sumMlbJESUpPdfVal", ROOT.AddressOf(data, "sumMlbJESUpPdfVal"), "sumMlbJESUpPdfVal/F")
        newTree.Branch("sumMlbJESDownPdfVal", ROOT.AddressOf(data, "sumMlbJESDownPdfVal"), "sumMlbJESDownPdfVal/F")
        newTree.Branch("deltaPhiPdfVal", ROOT.AddressOf(data, "deltaPhiPdfVal"), "deltaPhiPdfVal/F")
        newTree.Branch("metPdfVal", ROOT.AddressOf(data, "metPdfVal"), "metPdfVal/F")
        newTree.Branch("metJESUpPdfVal", ROOT.AddressOf(data, "metJESUpPdfVal"), "metJESUpPdfVal/F")
        newTree.Branch("metJESDownPdfVal", ROOT.AddressOf(data, "metJESDownPdfVal"), "metJESDownPdfVal/F")
        newTree.Branch("nLL", ROOT.AddressOf(data, "nLL"), "nLL/F")
        newTree.Branch("nLLJESUp", ROOT.AddressOf(data, "nLLJESUp"), "nLLJESUp/F")
        newTree.Branch("nLLJESDown", ROOT.AddressOf(data, "nLLJESDown"), "nLLJESDown/F")
        newTree.Branch("weight", ROOT.AddressOf(data, "weight"), "weight/F")
        newTree.Branch("prefireWeight", ROOT.AddressOf(data, "prefireWeight"), "prefireWeight/F")
        #~ newTree.Branch("weightUp", ROOT.AddressOf(data, "weightUp"), "weightUp/F")
        #~ newTree.Branch("weightDown", ROOT.AddressOf(data, "weightDown"), "weightDown/F")
        newTree.Branch("genWeight", ROOT.AddressOf(data, "genWeight"), "genWeight/F")
        newTree.Branch("bTagWeightLoose", ROOT.AddressOf(data, "bTagWeightLoose"), "bTagWeightLoose/F")
        newTree.Branch("bTagWeight", ROOT.AddressOf(data, "bTagWeight"), "bTagWeight/F")
        newTree.Branch("genPtTop1", ROOT.AddressOf(data, "genPtTop1"), "genPtTop1/F")
        newTree.Branch("genPtTop2", ROOT.AddressOf(data, "genPtTop2"), "genPtTop2/F")
        newTree.Branch("leptonFullSimScaleFactor1", ROOT.AddressOf(data, "leptonFullSimScaleFactor1"), "leptonFullSimScaleFactor1/F")
        newTree.Branch("leptonFullSimScaleFactor2", ROOT.AddressOf(data, "leptonFullSimScaleFactor2"), "leptonFullSimScaleFactor2/F")
        newTree.Branch("leptonFullSimScaleFactorErr1", ROOT.AddressOf(data, "leptonFullSimScaleFactorErr1"), "leptonFullSimScaleFactorErr1/F")
        newTree.Branch("leptonFullSimScaleFactorErr2", ROOT.AddressOf(data, "leptonFullSimScaleFactorErr2"), "leptonFullSimScaleFactorErr2/F")
        newTree.Branch("kfactor_zzto4l_pt", ROOT.AddressOf(data, "kfactor_zzto4l_pt"), "kfactor_zzto4l_pt/F")
        newTree.Branch("kfactor_zzto2l_pt", ROOT.AddressOf(data, "kfactor_zzto2l_pt"), "kfactor_zzto2l_pt/F")
        newTree.Branch("kfactor_zzto4l_mass", ROOT.AddressOf(data, "kfactor_zzto4l_mass"), "kfactor_zzto4l_mass/F")
        newTree.Branch("kfactor_zzto2l_mass", ROOT.AddressOf(data, "kfactor_zzto2l_mass"), "kfactor_zzto2l_mass/F")
        newTree.Branch("motherPdgId1", ROOT.AddressOf(data, "motherPdgId1"), "motherPdgId1/I")
        newTree.Branch("motherPdgId2", ROOT.AddressOf(data, "motherPdgId2"), "motherPdgId2/I")
        newTree.Branch("grandMotherPdgId1", ROOT.AddressOf(data, "grandMotherPdgId1"), "grandMotherPdgId1/I")
        newTree.Branch("grandMotherPdgId2", ROOT.AddressOf(data, "grandMotherPdgId2"), "grandMotherPdgId2/I")
        
        if tree == None or tree.GetEntries() == 0:
                return newTree
        
        from array import array
        weight = array('f', [0])
        tree.SetBranchAddress('weight', weight)
        genWeight = array('f', [0])
        tree.SetBranchAddress('genWeight', genWeight)
        prefireWeight = array('f', [0])
        tree.SetBranchAddress('prefireWeight', prefireWeight)
        bTagWeight = array('f', [0])
        tree.SetBranchAddress('bTagWeight', bTagWeight)
        if isMC:
                bTagWeightLoose = array('f', [0])
                tree.SetBranchAddress('bTagWeightLoose', bTagWeightLoose)
                kfactor_zzto4l_pt = array('f', [0])
                tree.SetBranchAddress('kfactor_zzto4l_pt', kfactor_zzto4l_pt)
                kfactor_zzto2l_pt = array('f', [0])
                tree.SetBranchAddress('kfactor_zzto2l_pt', kfactor_zzto2l_pt)
                kfactor_zzto4l_mass = array('f', [0])
                tree.SetBranchAddress('kfactor_zzto4l_mass', kfactor_zzto4l_mass)
                kfactor_zzto2l_mass = array('f', [0])
                tree.SetBranchAddress('kfactor_zzto2l_mass', kfactor_zzto2l_mass)
        genPtTop1 = array('f', [0])
        tree.SetBranchAddress('genPtTop1', genPtTop1)
        genPtTop2 = array('f', [0])
        tree.SetBranchAddress('genPtTop2', genPtTop2)
        leptonFullSimScaleFactor1 = array('f', [0])
        tree.SetBranchAddress('leptonFullSimScaleFactor1', leptonFullSimScaleFactor1)
        leptonFullSimScaleFactor2 = array('f', [0])
        tree.SetBranchAddress('leptonFullSimScaleFactor2', leptonFullSimScaleFactor2)
        leptonFullSimScaleFactorErr1 = array('f', [0])
        tree.SetBranchAddress('leptonFullSimScaleFactorErr1', leptonFullSimScaleFactorErr1)
        leptonFullSimScaleFactorErr2 = array('f', [0])
        tree.SetBranchAddress('leptonFullSimScaleFactorErr2', leptonFullSimScaleFactorErr2)
        motherPdgId1 = array('i', [0])
        tree.SetBranchAddress('motherPdgId1', motherPdgId1)
        motherPdgId2 = array('i', [0])
        tree.SetBranchAddress('motherPdgId2', motherPdgId2)
        grandMotherPdgId1 = array('i', [0])
        tree.SetBranchAddress('grandMotherPdgId1', grandMotherPdgId1)
        grandMotherPdgId2 = array('i', [0])
        tree.SetBranchAddress('grandMotherPdgId2', grandMotherPdgId2)
        mll = array('f', [0])
        tree.SetBranchAddress('mll', mll)
        chargeProduct = array('f', [0])
        tree.SetBranchAddress('chargeProduct', chargeProduct)
        nJets = array('i', [0])
        tree.SetBranchAddress('nJets', nJets)
        nShiftedJetsJESUp = array('i', [0])
        tree.SetBranchAddress('nShiftedJetsJESUp', nShiftedJetsJESUp)
        nShiftedJetsJESDown = array('i', [0])
        tree.SetBranchAddress('nShiftedJetsJESDown', nShiftedJetsJESDown)
        nBJets = array('i', [0])
        tree.SetBranchAddress('nBJets', nBJets)
        nBadMuonJets = array('i', [0])
        tree.SetBranchAddress('nBadMuonJets', nBadMuonJets)
        nLooseLeptons = array('i', [0])
        tree.SetBranchAddress('nLooseLeptons', nLooseLeptons)
        nVertices = array('i', [0])
        tree.SetBranchAddress('nVertices', nVertices)
        runNr = array('i', [0])
        tree.SetBranchAddress('runNr', runNr)
        lumiSec = array('i', [0])
        tree.SetBranchAddress('lumiSec', lumiSec)
        eventNr = array('l', [0])
        tree.SetBranchAddress('eventNr', eventNr)
        ht = array('f', [0])
        tree.SetBranchAddress('ht', ht)
        met = array('f', [0])
        tree.SetBranchAddress('met', met)
        metJESUp = array('f', [0])
        tree.SetBranchAddress('metJESUp', metJESUp)
        metJESDown = array('f', [0])
        tree.SetBranchAddress('metJESDown', metJESDown)
        caloMet = array('f', [0])
        tree.SetBranchAddress('caloMet', caloMet)
        genMet = array('f', [0])
        tree.SetBranchAddress('genMet', genMet)
        pt1 = array('f', [0])
        tree.SetBranchAddress('pt1', pt1)
        pt2 = array('f', [0])
        tree.SetBranchAddress('pt2', pt2)
        eta1 = array('f', [0])
        tree.SetBranchAddress('eta1', eta1)
        eta2 = array('f', [0])
        tree.SetBranchAddress('eta2', eta2)
        deltaPhiJetMet1 = array('f', [0])
        tree.SetBranchAddress('deltaPhiJetMet1', deltaPhiJetMet1)
        deltaPhiJetMet2 = array('f', [0])
        tree.SetBranchAddress('deltaPhiJetMet2', deltaPhiJetMet2)
        MT2 = array('f', [0])
        tree.SetBranchAddress('MT2', MT2)
        pt = array('f', [0])
        tree.SetBranchAddress('pt', pt)
        sumMlb = array('f', [0])
        tree.SetBranchAddress('sumMlb', sumMlb)
        sumMlbJESUp = array('f', [0])
        tree.SetBranchAddress('sumMlbJESUp', sumMlbJESUp)
        sumMlbJESDown = array('f', [0])
        tree.SetBranchAddress('sumMlbJESDown', sumMlbJESDown)
        deltaPhi = array('f', [0])
        tree.SetBranchAddress('deltaPhi', deltaPhi)
        deltaR = array('f', [0])
        tree.SetBranchAddress('deltaR', deltaR)

        
        # Initialize RooFit stuff that is slow
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
        
        lepsZPt_Edge = w.var("lepsZPt_Edge")
        sum_mlb_Edge =  w.var("sum_mlb_Edge")
        lepsDPhi_Edge =  w.var("lepsDPhi_Edge")
        met_Edge =  w.var("met_Edge")
        
        obsPt = ROOT.RooArgSet(lepsZPt_Edge)
        obsMlb = ROOT.RooArgSet(sum_mlb_Edge)
        obsDeltaPhi = ROOT.RooArgSet(lepsDPhi_Edge)
        obsMet = ROOT.RooArgSet(met_Edge)
        
        i = 0
        # Fill tree
        while tree.GetEntry(i):
                i+=1
                data.weight = weight[0]
                data.prefireWeight = prefireWeight[0]
                data.genWeight = genWeight[0]
                if isMC:
                        data.bTagWeight = bTagWeight[0]
                        data.bTagWeightLoose = bTagWeightLoose[0]
                else:
                        data.bTagWeight = 1.
                        data.bTagWeightLoose = 1.
                data.genPtTop1 = genPtTop1[0]
                data.genPtTop2 = genPtTop2[0]
                if not isMC:
                        data.leptonFullSimScaleFactor1 = 1.
                        data.leptonFullSimScaleFactor2 = 1.
                        data.leptonFullSimScaleFactorErr1 = 0.
                        data.leptonFullSimScaleFactorErr2 = 0.
                        data.kfactor_zzto4l_pt = 1.
                        data.kfactor_zzto2l_pt = 1.
                        data.kfactor_zzto4l_mass = 1.
                        data.kfactor_zzto2l_mass = 1.
                else:                                        
                        data.leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1[0]
                        data.leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2[0]
                        data.leptonFullSimScaleFactorErr1 = leptonFullSimScaleFactorErr1[0]
                        data.leptonFullSimScaleFactorErr2 = leptonFullSimScaleFactorErr2[0]
                        #data.leptonFullSimScaleFactorErr1 =  0
                        #data.leptonFullSimScaleFactorErr2 =  0
                        data.kfactor_zzto4l_pt = kfactor_zzto4l_pt[0]
                        data.kfactor_zzto2l_pt = kfactor_zzto2l_pt[0]
                        data.kfactor_zzto4l_mass = kfactor_zzto4l_mass[0]
                        data.kfactor_zzto2l_mass = kfactor_zzto2l_mass[0]
                        
                data.motherPdgId1 = motherPdgId1[0]
                data.motherPdgId2 = motherPdgId2[0]
                data.grandMotherPdgId1 = grandMotherPdgId1[0]
                data.grandMotherPdgId2 = grandMotherPdgId2[0]
                data.mll = mll[0]
                data.chargeProduct = chargeProduct[0]
                data.nJets = nJets[0]
                data.nShiftedJetsJESUp = nShiftedJetsJESUp[0]
                data.nShiftedJetsJESDown = nShiftedJetsJESDown[0]
                data.nBJets = nBJets[0]
                data.nBadMuonJets = nBadMuonJets[0]
                data.nLooseLeptons = nLooseLeptons[0]
                #data.nIsoTracks = tree.nIsoTracks
                data.nVertices = nVertices[0]
                data.runNr = runNr[0]
                data.lumiSec = lumiSec[0]
                data.eventNr = eventNr[0]
                data.vetoHEM = 1
                data.affectedHEM = 0
                if not isMC and year == 2018:
                        # hem treatment in data
                        # veto electrons and jets in affected region
                        data.affectedHEM = tree.affectedHEM
                        data.vetoHEM = tree.vetoHEM
                        #if runNr[0] >= 319077:
                                #data.affectedHEM = 1
                                #if abs(tree.pdgId1) == 11:
                                        #if tree.lepton1.Eta() > -4.7 and tree.lepton1.Eta() < -1.4 and tree.lepton1.Phi() > -1.6 and tree.lepton1.Phi() < -0.8:
                                               #data.vetoHEM = 0
                                #if abs(tree.pdgId2) == 11:
                                        #if tree.lepton2.Eta() > -4.7 and tree.lepton2.Eta() < -1.4 and tree.lepton2.Phi() > -1.6 and tree.lepton2.Phi() < -0.8:
                                               #data.vetoHEM = 0
                                #if tree.jet1.Eta() > -4.7 and tree.jet1.Eta() < -1.4 and tree.jet1.Phi() > -1.6 and tree.jet1.Phi() < -0.8:
                                        #data.vetoHEM = 0
                                #if tree.jet2.Eta() > -4.7 and tree.jet2.Eta() < -1.4 and tree.jet2.Phi() > -1.6 and tree.jet2.Phi() < -0.8:
                                        #data.vetoHEM = 0
                                
                        
                        pass
                elif year == 2018:
                        # hem treatment in MC
                        data.affectedHEM = tree.affectedHEM
                        data.vetoHEM = tree.vetoHEM
                        # veto fraction of events: 38.58 of 58.82fb^-1 lumi affected ~= 1286/1961
                        #fracNum = 1286
                        #fracDen = 1961
                        #if eventNr[0] % fracDen < fracNum:
                                #data.affectedHEM = 1
                                #if abs(tree.pdgId1) == 11:
                                       #if tree.lepton1.Eta() > -4.7 and tree.lepton1.Eta() < -1.4 and tree.lepton1.Phi() > -1.6 and tree.lepton1.Phi() < -0.8:
                                               #data.vetoHEM = 0
                                #if abs(tree.pdgId2) == 11:
                                       #if tree.lepton2.Eta() > -4.7 and tree.lepton2.Eta() < -1.4 and tree.lepton2.Phi() > -1.6 and tree.lepton2.Phi() < -0.8:
                                               #data.vetoHEM = 0
                                #if tree.jet1.Eta() > -4.7 and tree.jet1.Eta() < -1.4 and tree.jet1.Phi() > -1.6 and tree.jet1.Phi() < -0.8:
                                        #data.vetoHEM = 0
                                #if tree.jet2.Eta() > -4.7 and tree.jet2.Eta() < -1.4 and tree.jet2.Phi() > -1.6 and tree.jet2.Phi() < -0.8:
                                        #data.vetoHEM = 0
                data.ht = ht[0]
                data.met = met[0]
                data.metJESUp = metJESUp[0]
                data.metJESDown = metJESDown[0]
                data.caloMet = caloMet[0]
                data.genMet = genMet[0]
                data.pt1 = pt1[0]
                data.pt2 = pt2[0]
                data.eta1 = eta1[0]
                data.eta2 = eta2[0]
                data.deltaPhiJetMet1 = deltaPhiJetMet1[0]
                data.deltaPhiJetMet2 = deltaPhiJetMet2[0]
                data.MT2 = MT2[0]
                data.pt = pt[0]
                data.sumMlb = sumMlb[0]
                data.sumMlbJESUp = sumMlbJESUp[0]
                data.sumMlbJESDown = sumMlbJESDown[0]
                data.deltaPhi = abs(deltaPhi[0])
                data.deltaR = deltaR[0]
                
                lepsZPt_Edge.setVal(pt[0])
                sum_mlb_Edge.setVal(sumMlb[0])
                lepsDPhi_Edge.setVal(abs(deltaPhi[0]))
                met_Edge.setVal(met[0])
                

                data.ptPdfVal = ptPdf.getVal(obsPt)
                data.sumMlbPdfVal = sumMlbPdf.getVal(obsMlb)
                data.deltaPhiPdfVal = deltaPhiPdf.getVal(obsDeltaPhi)
                data.metPdfVal = metPdf.getVal(obsMet)
                
                met_Edge.setVal(metJESUp[0])
                sum_mlb_Edge.setVal(sumMlbJESUp[0])  
                data.metJESUpPdfVal = metPdf.getVal(obsMet)
                data.sumMlbJESUpPdfVal = sumMlbPdf.getVal(obsMlb)       
                met_Edge.setVal(metJESDown[0])
                sum_mlb_Edge.setVal(sumMlbJESDown[0])
                data.metJESDownPdfVal = metPdf.getVal(obsMet)   
                data.sumMlbJESDownPdfVal = sumMlbPdf.getVal(obsMlb)     
                
                data.nLL = - TMath.Log(data.ptPdfVal*data.sumMlbPdfVal*data.metPdfVal*data.deltaPhiPdfVal)
                data.nLLJESUp = - TMath.Log(data.ptPdfVal*data.sumMlbJESUpPdfVal*data.metJESUpPdfVal*data.deltaPhiPdfVal)
                data.nLLJESDown = - TMath.Log(data.ptPdfVal*data.sumMlbJESDownPdfVal*data.metJESDownPdfVal*data.deltaPhiPdfVal)
                newTree.Fill()
        
        pdfSourceFile.Close()
        w.Delete()
        return newTree

                                        

def main():
        
        
        

        parser = argparse.ArgumentParser(description='edge fitter reloaded.')
        
        parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False,
                                                  help="Verbose mode.")
        parser.add_argument("-S", "--sample", action="store", dest="sample", default="all",
                                                  help="choose a sample to convert, default is all.")
        parser.add_argument("-s", "--selection", dest = "selection" , action="store", default="SignalInclusive",
                                                  help="selection which to apply.")     
        parser.add_argument("-y", "--year", dest = "year" , action="store", default=2016, type=int,
                                                  help="Trees from which year to use.")     
        args = parser.parse_args()      
        
        cut = " chargeProduct < 0 && (pt1 > 25 || pt2 > 25)   && deltaR > 0.1 && p4.M() > 20"
        # && ((abs(eta1) < 1.4 || abs(eta1) > 1.6) && (abs(eta2) < 1.4 || abs(eta2) > 1.6)) && abs(eta1)<2.4  && abs(eta2) < 2.4 
        
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
                
        path = locations[args.year].dataSetPath
        outpath = locations[args.year].dataSetPathNLL
        if path[-1] == "/":
                swversion = path.split("/")[-2]
        else:
                swversion = path.split("/")[-1]
        cutsversion = "cuts%d"%(args.year)
        inputFiles = {
                                        "pdfs": {       2016:"workspace_NLLpdfs_sw2016v1019_SignalInclusive_Run2016_36fb.root",
                                                        2017:"workspace_NLLpdfs_sw2017v1005_SignalInclusive_Run2017_42fb.root",
                                                        2018:"workspace_NLLpdfs_sw2018v1001_SignalInclusive_Run2018_60fb.root",
                                                },
                                        }
        
        
        # init ROOT
        gROOT.Reset()
        
        ensurePathExists(outpath)
        
        samples = getFilePathsAndSampleNames(path)
        
        treeEEraw = readTrees(path, "EE")
        treeEMuraw = readTrees(path, "EMu")
        treeMuMuraw = readTrees(path, "MuMu")
        
        ### Only use events where the MET filters are ok
        baseCut = cut+ " && metFilterSummary > 0"
        
        cut = baseCut
        
        cutEE = cut + "&& triggerSummary > 0"
        cutEM = cut + "&& triggerSummary > 0"
        cutMM = cut + "&& triggerSummary > 0"
        
        for sample in samples:
                if sample == "WWTo4Q_Powheg_Spring16_25ns" or sample == "WWTo4Q_Powheg_Summer16_25ns" or sample == "WWTo4Q_Powheg_Fall17_25ns":  ### No events remain
                        continue
                if sample == "MergedData":
                        isMC = False
                else:
                        isMC = True
                if args.sample != "all" and args.sample not in sample:
                        continue
                
                inFile = TFile(samples[sample], "read")
                denominator = inFile.FindObjectAny("analysis paths").GetBinContent(1)
                histo = TH1F("analysis paths", "analysis paths", 2, 0 ,1)
                #histo.SetDirectory(0)
                histo.GetXaxis().SetBinLabel(1, "None")
                histo.SetBinContent(1, denominator)
                
                
                print sample
                f1 = TFile("%s/%s.processed.%s.root"%(outpath,swversion,sample),"RECREATE")
                folder = f1.mkdir("%sDileptonFinalTrees"%cutsversion)
                

                treeEEConverted = convertDileptonTree(treeEEraw[sample].CopyTree(cutEE),inputFiles,isMC,args.year)
                treeEEraw[sample] = None
                
                treeEMuConverted = convertDileptonTree(treeEMuraw[sample].CopyTree(cutEM),inputFiles,isMC,args.year)
                treeEMuraw[sample] = None

                treeMuMuConverted = convertDileptonTree(treeMuMuraw[sample].CopyTree(cutMM),inputFiles,isMC,args.year)
                treeMuMuraw[sample] = None
                
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
                inFile.Close()
                treeEEConverted.IsA().Destructor(treeEEConverted)  
                treeEMuConverted.IsA().Destructor(treeEMuConverted)        
                treeMuMuConverted.IsA().Destructor(treeMuMuConverted)  
                

                
                        
                        
        
                        

main()
