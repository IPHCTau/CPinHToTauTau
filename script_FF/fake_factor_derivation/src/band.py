#! /usr/bin/env python3
# Example of fit with pol3 function F(x) = p0 + p1*x + p2*x^2 + p3*x^3
# showing analytic expression for uncertainty bands 

import ROOT
from array import array

# fitted function (pol3 : p0+p1*x+p2*x^2+p3*x^3
def FitFunc(x,par):
    f = 0
    xx = 1
    for i in range(0,4):
        f += par[i]*xx
        xx *= x[0]
    return f
    
# uncertainty band
def Band(x,par):
    W = 0
    xx = 1
    for i in range(0,7):
        W += par[i]*xx
        xx *= x[0]
    return ROOT.TMath.Sqrt(W)

# upper uncertainty 
def FitFuncPlus(x,par):
    parBand = []
    # shift index parameter by 4
    for i in range(0,7):
        parBand.append(par[i+4])
    return FitFunc(x,par)+Band(x,parBand)

# lower uncertainty 
def FitFuncMinus(x,par):
    parBand = []
    # shift index parameter by 2
    for i in range(0,7):
        parBand.append(par[i+4])
    return FitFunc(x,par)-Band(x,parBand)


def FitFuncMinus_str():
    
    # Construct the polynomial part (F(x) = p0 + p1*x + p2*x^2 + p3*x^3)
    poly_str = f"[0] + [1]*x + [2]*x**2 + [3]*x**3"
    
    # Construct the uncertainty band part (Band(x) = sqrt(par4 + par5*x + ... + par10*x^6))
    band_str = f"sqrt([4] + [5]*x + [6]*x**2 + [7]*x**3 + [8]*x**4 + [9]*x**5 + [10]*x**6)"
    
    # Return the full formula as a string
    return f"({poly_str}) - {band_str}"

def FitFuncPlus_str():
    poly_str = f"[0] + [1]*x + [2]*x**2 + [3]*x**3"
    band_str = f"sqrt([4] + [5]*x + [6]*x**2 + [7]*x**3 + [8]*x**4 + [9]*x**5 + [10]*x**6)"
    return f"({poly_str}) + {band_str}"


    
############
#   MAIN   #
############
def fit_upper_lower(fit, fit_result,h_uncert):

    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)

    

    # fitted function (p0+p1*x+p2*x^2+p3*x^3)
    # 4 parameters
    fitFunc = fit
    # some initial parameters ->
    fitFunc.SetLineColor(ROOT.kBlue)


    FitFuncPlus_str_ = FitFuncPlus_str()
    FitFuncMinus_str_ = FitFuncMinus_str()
    # upper uncertainty
    # 4 + 7 = 11 parameters
    fitFuncPlus = ROOT.TF1('FitFuncPlus',FitFuncPlus_str_,35,200.,11)
    fitFuncPlus.SetLineColor(ROOT.kRed)

    # lower uncertainty
    fitFuncMinus = ROOT.TF1('FitFuncMinus',FitFuncMinus_str_,35.,200.,11)
    fitFuncMinus.SetLineColor(ROOT.kRed)
    
    canv = ROOT.TCanvas('canv','canv',600,600)
    fitResult = fit_result
    # getting covariance matrix ->
    # cov(i,j) = err(i) * err(j) * rho(i,j)
    # here err(i) is the uncertainty in the fitted parameter "i"
    # and rho(i,j) is the correlation matrix of the fitted parameters
    cov = fitResult.GetCovarianceMatrix()
    # Parameters of the uncertainty band functions
    fitFuncPlus.SetParameter(0,fitFunc.GetParameter(0))
    fitFuncPlus.SetParameter(1,fitFunc.GetParameter(1))
    fitFuncPlus.SetParameter(2,fitFunc.GetParameter(2))
    fitFuncPlus.SetParameter(3,fitFunc.GetParameter(3))
    fitFuncPlus.SetParameter(4,cov(0,0))
    fitFuncPlus.SetParameter(5,2*cov(0,1))
    fitFuncPlus.SetParameter(6,2*cov(0,2)+cov(1,1))
    fitFuncPlus.SetParameter(7,2*cov(0,3)+2*cov(1,2))
    fitFuncPlus.SetParameter(8,2*cov(1,3)+cov(2,2))
    fitFuncPlus.SetParameter(9,2*cov(2,3))
    fitFuncPlus.SetParameter(10,cov(3,3))
    for i in range(0,11):
        fitFuncMinus.SetParameter(i,fitFuncPlus.GetParameter(i))


    # ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hfit,0.68)
    # hfit.Draw('e2') # from TVirtualFitter
    h_uncert.Draw('e2') # from TVirtualFitter
    # cyan transparent band
    h_uncert.SetFillColor(ROOT.kCyan)
    h_uncert.SetFillStyle(3001)
    h_uncert.SetMarkerSize(0)
    h_uncert.SetLineWidth(0)

    fit.Draw('lsame') # fitted function
    fitFuncPlus.Draw('lsame') # analytic function for upper band
    fitFuncMinus.Draw('lsame') # analytic function for lower band
    # hist.Draw('e1same') # fitted measurement points

    print("fitFuncPlus: ", fitFuncPlus)
    print("fitFuncMinus: ", fitFuncMinus)

    print("fitFuncPlus formula: ", fitFuncPlus.GetExpFormula('P'))
    print("fitFuncMinus formula: ", fitFuncMinus.GetExpFormula('P'))
    
    canv.RedrawAxis()
    canv.Update()
    canv.Print('example_pol3.png')
    canv.SaveAs('example_pol3.root')



    return fitFuncPlus, fitFuncMinus