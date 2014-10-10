
from ROOT import gRandom, TH1, TH1D, gROOT, cout

gROOT.LoadMacro( "/home/skluth/Downloads/RooUnfold/libRooUnfold.so" )

from ROOT import RooUnfoldResponse
from ROOT import RooUnfoldBasisSplines

from ROOT import TProfile, TCanvas, TF1, TMatrixD, TH2D, TVectorD

TH1.SetDefaultSumw2( False )



def accept( x ):
    p= 1.0 - 0.5*(x-1.0)**2
    r= gRandom.Rndm()
    result= True
    if r > p: 
        result= False
    return result

def transform( x ):
    return x - 0.2*x**2/4.0

def smear( x, mean=0.0, sigma=0.1 ):
    xsmear= gRandom.Gaus( mean, sigma )
    return x + xsmear

def trueFun( xx, par ):
    xxx= xx[0]
    b= [ 1.0, 10.0, 5.0 ]
    x= [ 0.4, 0.8, 1.5 ]
    g= [ 2.0, 0.2, 0.2 ]
    f= 0.0
    for k in range( 3 ):
        f+= b[k]*g[k]**2/(( xxx - x[k] )**2+g[k]**2)
    return f*par[0]


def run():

    global tf, tf2, htrue, htrue2, hmeas, hreco, canv

    tf= TF1( "tf", trueFun, 0.0, 2.0, 1 )
    tf.SetParameter( 0, 1.0 )

    tf2= TF1( "tf2", trueFun, 0.0, 2.0, 1 )

    tfint= tf.Integral( 0.0, 2.0 )
    nevnt= 5000

    htrue= TH1D( "htrue", "True distribution", 40, 0.0, 2.0 )
    htrue2= TH1D( "htrue2", "True distribution", 12, 0.0, 2.0 )
    hmeas= TH1D( "hmeas", "Measured distribution", 40, 0.0, 2.0 )
    res= RooUnfoldResponse( 40, 0.0, 2.0, 12, 0.0, 2.0 )

    for i in xrange( nevnt ):
        f= tf.GetRandom()
        htrue.Fill( f )
        htrue2.Fill( f )
        if accept( f ):
            fsmear= smear( transform( f ) )
            hmeas.Fill( fsmear )
            res.Fill( fsmear, f )
        else:
            res.Miss( f )

    canv= TCanvas( "canv", "Blobel unfolding example", 600, 800 )

    unfold= RooUnfoldBasisSplines( res, hmeas, 1, 0.0, 0, 2 )
    unfold.SetVerbose( 2 )
    unfold.PrintTable( cout, htrue2, 2 )
    hreco= unfold.Hreco( 2 )
    
    scale= float(nevnt)/tfint*(2.0/40.0)
    tf.SetParameter( 0, scale )

    canv.Divide( 1, 3 )

    canv.cd( 1 )
    htrue.Draw( "e" )
    tf.Draw( "same" )

    canv.cd( 2 )
    hmeas.Draw( "e" )
    tf.Draw( "same" )
    
    scale= float(nevnt)/tfint*(2.0/12.0)
    tf2.SetParameter( 0, scale )
 
    canv.cd( 3 )
    hreco.SetMaximum( 1000.0 )
    hreco.Draw( "e" )
    htrue2.Draw( "same" )
    tf2.Draw( "same" )
    

    return
