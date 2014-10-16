#!/usr/bin/env python

from ROOT import gRandom, TH1, TH1D, gROOT, cout, TMath, gStyle

gROOT.LoadMacro( "$HOME/Downloads/RooUnfold/libRooUnfold.so" )

from ROOT import RooUnfoldResponse
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import RooUnfoldTUnfold
from ROOT import RooUnfoldInvert
from ROOT import RooUnfoldBasisSplines

from ROOT import TProfile, TCanvas, TF1, TMatrixD, TH2D, TVectorD

TH1.SetDefaultSumw2( False )


#  Gaussian smearing, systematic translation, and variable inefficiency:

class Measurement:

    def __init__( self, mean, sigma, leff ):
        self.mean= mean
        self.sigma= sigma
        self.leff= leff
        return

    def getMean( self ):
        return self.mean
    def getSigma( self ):
        return self.sigma
    def getLeff( self ):
        return self.leff

    def efficiency( self, x ):
        return 0.3 + (1.0-0.3)/10.0*x

    def accept( self, x ):
        result= True
        if self.leff:
            p= self.efficiency( x )
            r= gRandom.Rndm()
            if r > p:
                result= False
        return result

    def smear( self, xt ):
        xsmear= gRandom.Gaus( self.mean, self.sigma )
        return xt+xsmear

    def transform( self, x ):
        return x

    def generate( self, fun ):
        xt= fun.GetRandom()
        x= None
        if self.accept( xt ):
            x= self.smear( self.transform( xt ) )
        return xt, x

class BlobelMeasurement( Measurement ):
    def __init__( self, mean, sigma, leff ):
        Measurement.__init__( self, mean, sigma, leff )
        return
    def efficiency( self, x ):
        x= x/5.0
        return 1.0 - 0.5*(x-1.0)**2
    def transform( self, x ):
        x= x/5.0
        return ( x - 0.2*x**2/4.0 )*5.0

def createMeasurement( mean, sigma, leff, opt="" ):
    if "blobel" in opt:
        return BlobelMeasurement(  mean, sigma, leff )
    else:
        return Measurement(  mean, sigma, leff )

class Trainer:

    def __init__( self, bininfo, opttfun, optfun ):
        self.bininfo= bininfo.create( optfun )
        self.trainfun= createFunction( self.bininfo, opttfun, "train_"+opttfun )
        self.opttfun= opttfun
        return

    # Training: create response object:
    def train( self, measurement, loufl=False ):
        print "============================ TRAIN ============================="
        txt= "Smear mu, s.d.: " + str(measurement.getMean())
        txt+= ", " + str(measurement.getSigma()) + ", eff.: "
        txt+= str(measurement.getLeff())
        txt+= ", o/u-flow: " + str(loufl) + ", function: " + self.opttfun
        print txt
        # Create response matrix object:
        response= RooUnfoldResponse( self.bininfo["mbins"],
                                     self.bininfo["mlo"],
                                     self.bininfo["mhi"],
                                     self.bininfo["tbins"]*self.bininfo["nrebin"],
                                     self.bininfo["tlo"],
                                     self.bininfo["thi"] )
        response.UseOverflow( loufl )
        for i in xrange( 100000 ):
            xt, x= measurement.generate( self.trainfun )
            if x != None:
                response.Fill( x, xt )
            else:
                response.Miss( xt )
        return response

class Tester:

    def __init__( self, bininfo, optfun ):
        self.bininfo= bininfo.create( optfun )
        self.testfun= createFunction( self.bininfo, optfun, "test_"+optfun )
        self.optfun= optfun
        return

    # Generate test distributions:
    def generateTest( self, measurement ):
        hTrue, hMeas= self.makeTestHistos()
        for i in xrange( 10000 ):
            xt, x= measurement.generate( self.testfun )
            hTrue.Fill( xt )
            if x != None:
                hMeas.Fill( x )
        return hTrue, hMeas

    def makeTestHistos( self ):
        hTrue= TH1D( "true", "Test Truth",
                     self.bininfo["tbins"]*self.bininfo["nrebin"],
                     self.bininfo["tlo"],
                     self.bininfo["thi"] )
        hMeas= TH1D( "meas", "Test Measured",
                     self.bininfo["mbins"],
                     self.bininfo["mlo"],
                     self.bininfo["mhi"] )
        return hTrue, hMeas

class UnfoldTester:

    def __init__( self, optunf, nrebin=1 ):
        optunfs= [ "Bayes", "SVD", "TUnfold", "Invert", "Reverse",
                   "BasisSplines" ]
        if not optunf in optunfs:
            txt= "UnfoldTester: Unfolding option " + optunf + " not recognised" 
            raise ValueError( txt )
        self.optunf= optunf
        self.nrebin= nrebin
        return

    # Create unfolder object:
    def unfolderFactory( self, response, hMeas, BasisSpline_m0=0 ):
        if "Bayes" in self.optunf:
            # # Bayes unfolding until chi^2 cut (max 10):
            # unfold= RooUnfoldBayes( response, hMeas, 10, False, True )
            # Bayes unfolding (default 4 iterations):
            unfold= RooUnfoldBayes( response, hMeas )
        elif "SVD" in self.optunf:
            # SVD unfolding with free regularisation:
            unfold= RooUnfoldSvd( response, hMeas )
        elif "TUnfold" in self.optunf:
            # TUnfold with fixed regularisation tau=0.002
            # unfold= RooUnfoldTUnfold( response, hMeas, 0.002 )
            unfold= RooUnfoldTUnfold( response, hMeas )
        elif "Invert" in self.optunf:
            unfold= RooUnfoldInvert( response, hMeas )
        elif "Reverse" in self.optunf:
            # Use equivalent Bayes with 1 iteration:
            unfold= RooUnfoldBayes( response, hMeas, 1 )
        elif "BasisSplines" in self.optunf:
            # BasisSplines with automatic tau determination
            unfold= RooUnfoldBasisSplines( response, hMeas, self.nrebin,
                                           0.0, BasisSpline_m0, 2 )
        return unfold

    # Run a test of the unfolding methods:
    def rununfoldtest( self, tester, measurement, response ):
        print "============================= TEST =========================="
        txt= "Method: " + self.optunf + ", smear mu, s.d.: "
        txt+= str(measurement.getMean()) + ", "
        txt+= str(measurement.getSigma())
        txt+= ", eff.: " + str(measurement.getLeff())
        txt+= ", function: " + tester.optfun
        print txt
        hTrue, hMeas= tester.generateTest( measurement )
        print "=========================== UNFOLD ========================"
        print "Unfolding method:", self.optunf
        unfold= self.unfolderFactory( response, hMeas )
        return unfold, hTrue, hMeas


# Create functions for training and testing:
def createFunction( bininfo, optfun, name ):
    if "exp" in optfun:
        fun= TF1( name, "exp(-x/3.0)",
                  bininfo["tlo"], bininfo["thi"] )
    elif "bw" in optfun:
        fun= TF1( name, "TMath::BreitWigner(x,4.0,1.0)",
                  bininfo["tlo"], bininfo["thi"] )
    elif "gaus" in optfun:
        fun= TF1( name, "TMath::Gaus(x,5.0,2.5)",
                  bininfo["tlo"], bininfo["thi"] )
    elif "blobel" in optfun:
        def blobelFun( xx ):
            # xxx= xx[0]
            xxx= xx[0]/5.0
            b= [ 1.0, 10.0, 5.0 ]
            x= [ 0.4, 0.8, 1.5 ]
            g= [ 2.0, 0.2, 0.2 ]
            f= 0.0
            for k in range( 3 ):
                f+= b[k]*g[k]**2/(( xxx - x[k] )**2+g[k]**2)
            return f
        fun= TF1( name, blobelFun, bininfo["tlo"], bininfo["thi"] )
    else:
        print "Function option", optfun, "not recognised"
        fun= None
    return fun

# Provide binning for different test functions:
class BinInfo:
    def __init__( self, nrebin=1 ):
        self.nrebin= nrebin
    def create( self, optfun ):
        if "exp" in optfun:
            bininfo= { "tbins": 10,
                       "tlo": 0.0,
                       "thi": 10.0,
                       "mbins": 20,
                       "mlo": -2.0,
                       "mhi": 10.0 }
        elif "gaus" in optfun:
            bininfo= { "tbins": 10,
                       "tlo": 0.0,
                       "thi": 10.0,
                       "mbins": 20,
                       "mlo": -2.0,
                       "mhi": 10.0 }
        elif "bw" in optfun:
            bininfo= { "tbins": 40,
                       "tlo": 0.0,
                       "thi": 10.0,
                       "mbins": 80,
                       "mlo": -2.0,
                       "mhi": 10.0 }
        elif "blobel" in optfun:
            bininfo= { "tbins": 12,
                       "tlo": 0.0,
                       # "thi": 2.0,
                       "thi": 10.0,
                       "mbins": 40,
                       "mlo": 0.0,
                       # "mhi": 2.0 }
                       "mhi": 10.0 }
        bininfo["nrebin"]= self.nrebin
        return bininfo

# Create text for histo titles:
def funtxts( opttfun, optfun, leff, loufl ):
    if "exp" in opttfun:
        funttxt= "Exp tau=3"
    elif "bw" in opttfun:
        funttxt= "B-W mu=4, s.d.=1"
    elif "gaus" in opttfun:
        funttxt= "Gaus mu=5, s.d.=2.5"
    elif "blobel" in opttfun:
        funttxt= "Blobel"
    if "exp" in optfun:
        funtxt= "Exp tau=3"
    elif "bw" in optfun:
        funtxt= "B-W mu=4, s.d.=1"
    elif "gaus" in optfun:
        funtxt= "Gaus mu=5, s.d.=2.5"
    elif "blobel" in optfun:
        funtxt= "Blobel"
    if leff:
        funtxt+= ", eff."
    if loufl:
        funtxt+= ", o/u-flow"
    return funttxt, funtxt

# Plot pull distributions:
def plotPulls( optunf="Bayes", ntest=10, leff=True, loufl=False,
               optfun="exp", opttfun="", gmean=-1.0, nrebin=4 ):

    if opttfun == "":
        opttfun= optfun

    if optfun == "blobel":
        gmean= 0.0

    funttxt, funtxt= funtxts( opttfun, optfun, leff, loufl )

    global histos, canv, canv2
    histos= []

    canv= TCanvas( "canv", "thruth vs reco pulls", 600, 800 )
    canv.Divide( 1, 3 )
    canv2= TCanvas( "canv2", "P(chi^2)", 600, 800 )
    canv2.Divide( 2, 3 )

    if optunf == "BasisSplines":
        bininfo= BinInfo( nrebin )
        unfoldtester= UnfoldTester( optunf, nrebin )
    else:
        bininfo= BinInfo()
        unfoldtester= UnfoldTester( optunf )
    trainer= Trainer( bininfo, opttfun, optfun )
    tester= Tester( bininfo, optfun )
    hbininfo= bininfo.create( optfun )
    dx= hbininfo["mhi"]
    for sigma, ipad in [ [ 0.01*dx, 1 ], [ 0.03*dx, 2 ], [ 0.1*dx, 3 ] ]:
        txt= optunf + ", smear mu, s.d.= " + str(gmean) + ", " + str(sigma) + ", train: " + funttxt + ", test: " + funtxt + ", " + str(ntest) + " tests"
        hPulls= TProfile( "pulls", txt,
                          hbininfo["tbins"], hbininfo["tlo"], hbininfo["thi"] )
        hPulls.SetErrorOption( "s" )
        hPulls.SetYTitle( "Thruth reco pull" )
        histos.append( hPulls )
        hChisq= TH1D( "chisq", "P(chi^2) rec " + txt, 10, 0.0, 1.0 )
        hChisqm= TH1D( "chisqm", "P(chi^2) mea " + txt, 10, 0.0, 1.0 )
        histos.append( hChisq )
        histos.append( hChisqm )
        measurement= createMeasurement( gmean, sigma, leff, optfun )        
        response= trainer.train( measurement, loufl=loufl )
        for itest in range( ntest ):
            print "Test", itest
            unfold, hTrue, hMeas= unfoldtester.rununfoldtest( tester,
                                                              measurement,
                                                              response )
            unfold.PrintTable( cout, hTrue, 2 )
            hReco= unfold.Hreco( 2 )
            nbin= hReco.GetNbinsX()
            if hbininfo["nrebin"] > 1:
                hTrue= hTrue.Rebin( nrebin )
            for ibin in range( nbin+1 ):
                truevalue= hTrue.GetBinContent( ibin )
                recvalue= hReco.GetBinContent( ibin )
                error= hReco.GetBinError( ibin )
                if error > 0.0:
                    pull= ( recvalue - truevalue )/error
                    hPulls.Fill( hReco.GetBinCenter( ibin ), pull )
            chisq= unfold.Chi2( hTrue, 2 )
            hChisq.Fill( TMath.Prob( chisq, hTrue.GetNbinsX() ) )
            chisqm= unfold.Chi2measured()
            pchisqm= TMath.Prob( chisqm, hMeas.GetNbinsX()-hReco.GetNbinsX() )
            print "Chisq measured=", chisqm, "P(chi^2)=", pchisqm
            hChisqm.Fill( pchisqm )

        canv.cd( ipad )
        gStyle.SetErrorX( 0 )
        hPulls.SetMinimum( -3.0 )
        hPulls.SetMaximum( 3.0 )
        hPulls.SetMarkerSize( 1.0 )
        hPulls.SetMarkerStyle( 20 )
        hPulls.SetStats( False )
        hPulls.Draw()
        canv2.cd( ipad*2-1 )
        hChisq.Draw()
        canv2.cd( ipad*2 )
        hChisqm.Draw()


    fname= "RooUnfoldTestPulls_" + optunf + "_" + opttfun + "_" + optfun
    if loufl:
        fname+= "_oufl"
    canv.Print( fname + ".pdf" )
    fname= "RooUnfoldTestChisq_" + optunf + "_" + opttfun + "_" + optfun
    if loufl:
        fname+= "_oufl"
    canv2.Print( fname + ".pdf" )
        
    return

# Plots with varying smearing resolution:
def featureSizePlots( optunf="Bayes", leff=True, optfun="exp", opttfun="",
                      loufl=False, gmean=-1.0, nrebin=4 ):

    if opttfun == "":
        opttfun= optfun

    if optfun == "blobel":
        gmean= 0.0

    funttxt, funtxt= funtxts( opttfun, optfun, leff, loufl )

    global hReco, hMeas, hTrue, hPulls, canv, histos
    histos= []
    canv= TCanvas( "canv", "feature size plots", 600, 800 )
    canv.Divide( 1, 3 )

    if optunf == "BasisSplines":
        bininfo= BinInfo( nrebin )
        unfoldtester= UnfoldTester( optunf, nrebin )
    else:
        bininfo= BinInfo()
        unfoldtester= UnfoldTester( optunf )
    trainer= Trainer( bininfo, opttfun, optfun )
    tester= Tester( bininfo, optfun )
    hbininfo= bininfo.create( optfun )
    dx= hbininfo["mhi"]
    for sigma, ipad in [ [ 0.01*dx, 1 ], [ 0.03*dx, 2 ], [ 0.1*dx, 3 ] ]:
        measurement= createMeasurement( gmean, sigma, leff, optfun )
        response= trainer.train( measurement, loufl=loufl )
        unfold, hTrue, hMeas= unfoldtester.rununfoldtest( tester,
                                                          measurement,
                                                          response )
        hReco= unfold.Hreco( 2 )

        hRecoMeasured= unfold.HrecoMeasured()

        chisqm= unfold.Chi2measured()
        pchisqm= TMath.Prob( chisqm, hMeas.GetNbinsX()-hReco.GetNbinsX() )
        print "Chisq measured=", chisqm, "P(chi^2)=", pchisqm

        if hbininfo["nrebin"] > 1:
            hTrue= hTrue.Rebin( nrebin )
        histos.append( [ hTrue, hMeas, hReco ] )
        canv.cd( ipad )
        gStyle.SetErrorX( 0 )
        hReco.SetTitle( optunf + ", smear mu, s.d.= " +
                        str(gmean) + ", " + str(sigma) +
                        ", train: " + funttxt + ", test: " + funtxt )
        hReco.SetYTitle( "Entries" )
        hReco.SetMinimum( 0.0 )
        hReco.SetMarkerStyle( 20 )
        hReco.SetMarkerSize( 0.75 )
        hReco.SetStats( False )
        hReco.Draw( "pe" )
        hRecoMeasured.Draw( "same" )
        hMeas.SetMarkerStyle( 20 )
        hMeas.SetMarkerSize( 0.5 )
        hMeas.Draw( "samepe" )
        hTrue.SetLineColor( 8 )
        hTrue.Draw( "same" )

    fname= "RooUnfoldTest_" + optunf + "_" + opttfun + "_" + optfun
    if loufl:
        fname+= "_oufl"
    canv.Print( fname + ".pdf" )

    return

# Make all plots of unfolding tests:
def doAllPlots( optunfs= "Bayes SVD TUnfold Invert Reverse BasisSplines",
                ntest=10 ):
    optunftokens= optunfs.split()
    for optunf in optunftokens:
        for optfun in [ "bw", "exp", "blobel" ]:
            for opttfun in [ "", "gaus" ]:
                # for loufl in [ False, True ]:
                #     featureSizePlots( optunf=optunf,
                #                       optfun=optfun, opttfun=opttfun,
                #                       loufl=loufl )
                #     plotPulls( optunf=optunf, ntest=ntest,
                #                optfun=optfun, opttfun=opttfun,
                #                loufl=loufl )
                featureSizePlots( optunf=optunf,
                                  optfun=optfun, opttfun=opttfun )
                plotPulls( optunf=optunf, ntest=ntest,
                           optfun=optfun, opttfun=opttfun )


    return



    

# Run as script:
if __name__ == '__main__':
    doAllPlots()

