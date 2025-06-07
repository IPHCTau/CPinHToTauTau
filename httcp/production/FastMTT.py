"""
  FastMTT : WiktorMat
  # https://github.com/WiktorMat/ClassicSVfit/blob/fastMTT_2024/python/FastMTT.py
  ### Main reference: https://github.com/SVfit/ClassicSVfit/blob/fastMTT_2024/src/FastMTT.cc
"""

import os

import law
import time
import functools
from typing import Optional
from columnflow.util import maybe_import

np = maybe_import("numpy")
#sc = maybe_import("scipy")
ak = maybe_import("awkward")
#mp = maybe_import("matplotlib")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
#physical_constants = sc.constants.physical_constants

logger = law.logger.get_logger(__name__)

ElectronMass = 0.51099895/1000  #MeV -> GeV 
MuonMass     = 105.6583755/1000 #MeV -> GeV 
ChargedPionMass = 139.5/1000 #MeV -> GeV 

#Invariant mass calculation
def InvariantMass(aP4):
    energy_squared = aP4[..., 3]**2
    momentum_squared = aP4[..., 0]**2 + aP4[..., 1]**2 + aP4[..., 2]**2
    return np.sqrt(energy_squared - momentum_squared)

class Likelihood:
    def __init__(self,
                 enable_MET = True,
                 enable_mass = True,
                 enable_BW = True,
                 enable_px = False,
                 enable_py = False,
                 enable_window=False,
                 constrain_window=[123.0, 127.0],
                 enable_gauss=False):

        #METinputs
        self.recoMET = np.array([0.0, 0.0, 0.0, 0.0])
        self.covMET = np.ones((2, 2))

        #setParameters
        self.coeff1 = 6
        self.coeff2 = 1/1.15

        #LeptonInputs
        self.leg1P4 = np.array([0.0, 0.0, 0.0, 0.0])
        self.leg2P4 = np.array([0.0, 0.0, 0.0, 0.0])

        #Visible mass of both leptons
        self.mvis = np.array([0.0, 0.0, 0.0, 0.0])

        #Invariant mass of each lepton
        self.mvisleg1 = np.array([0.0])
        self.mvisleg2 = np.array([0.0])

        self.mVisOverTauSquare1 = np.array([0.0])
        self.mVisOverTauSquare2 = np.array([0.0])
         
        self.mTau = 1776.86/1000 #MeV -> GeV
        
        self.leg1DecayType = np.array([0.0])
        self.leg2DecayType = np.array([0.0])
        self.leg1DecayMode = np.array([0.0])
        self.leg2DecayMode = np.array([0.0])


        #Enable/disable likelihood channel
        self.enable_MET = enable_MET
        self.enable_mass = enable_mass

        #These are experimental and not used by main code
        self.enable_BW = enable_BW
        self.enable_px = enable_px
        self.enable_py = enable_py
        self.enable_window = enable_window
        self.enable_gauss = enable_gauss

        self.window = constrain_window

        return

    def setLeptonInputs(self, aLeg1P4, aLeg2P4, aLeg1DecayType, aLeg2DecayType, aLeg1DecayMode, aLeg2DecayMode):
        
        self.leg1DecayType = aLeg1DecayType
        self.leg2DecayType = aLeg2DecayType
        self.leg1DecayMode = aLeg1DecayMode
        self.leg2DecayMode = aLeg2DecayMode

        self.leg1P4 = aLeg1P4
        self.leg2P4 = aLeg2P4
        
        #visible invariant mass
        #eq. (4)
        self.mvis = InvariantMass(self.leg1P4 + self.leg2P4)
        
        self.mvisleg1[(aLeg1DecayType==1) & (self.mvisleg1>1.5)] = 0.3
        self.mvisleg2[(aLeg2DecayType==1) & (self.mvisleg2>1.5)] = 0.3

        self.mVisOverTauSquare1 = (self.mvisleg1/self.mTau)**2
        self.mVisOverTauSquare2 = (self.mvisleg2/self.mTau)**2

    def massLikelihood(self, m: np.ndarray):
        mScaled = m*self.coeff2

        mask1 = (mScaled < self.mvis[:, np.newaxis])
        
        mVS2 = (self.mvis[:, np.newaxis]/mScaled)**2
        
        x1Min = np.minimum(1.0, self.mVisOverTauSquare1)
        x2Min = np.maximum(self.mVisOverTauSquare2[:, np.newaxis], mVS2)
        x2Max = np.minimum(1.0, mVS2/x1Min[:, np.newaxis])
        
        mask2 = (x2Min > x2Max)
        
        jacobiFactor = 2.0*self.mvis[:, np.newaxis]**2*mScaled**(-self.coeff1)
        x2IntegralTerm = np.log(x2Max/x2Min)

        value = 0.0
        value += x2IntegralTerm

        HadDecay1 = np.broadcast_to((self.leg1DecayType != 1)[:, np.newaxis], value.shape)
        value += HadDecay1 * mVS2 * (1 / x2Max - 1 / x2Min)

        HadDecay2 = np.broadcast_to((self.leg2DecayType != 1)[:, np.newaxis], value.shape)
        value += HadDecay2 * (mVS2*x2IntegralTerm - (x2Max - x2Min))

        value[mask1 | mask2] = 0.0

        value *= 1E9*jacobiFactor
        return value

    ###WORK IN PROGRESS###
    #This experimental component is still to be tested (and set to false by default)#
    #It will better constraint the likelihood function to Z0/H mass
    #(in order for better momenta estimation)

    def BreitWigner(self, invariant_mass):
        Higgs_mass = 125
        Higgs_gamma = Higgs_mass*0.01 #value set in original SVfit paper, however we will play with it yet
        Z_mass = 91.2
        Z_gamma = 2.5

        def normalization_constant(mass, gamma):
            x = mass*np.sqrt(mass**2+gamma**2)
            return 2*np.sqrt(2)*mass*gamma*x/np.pi/np.sqrt(mass**2 + x)

        H_denominator = (invariant_mass**2 - Higgs_mass**2)**2 + (Higgs_mass**2)*(Higgs_gamma**2)
        #Z_denominator = (invariant_mass**2 - Z_mass**2)**2 + (Z_mass**2)*(Z_gamma**2)
        #return normalization_constant(Z_mass, Z_gamma)/Z_denominator + normalization_constant(Higgs_mass, Higgs_gamma)/H_denominator
        return normalization_constant(Higgs_mass, Higgs_gamma)/H_denominator

    
    def Gauss(self, invariant_mass):
        Higgs_mass = 125
        Higgs_gamma = Higgs_mass*0.01 #value set in original SVfit paper, however we will play with it yet

        Higgs_gauss_factor = np.exp(-0.5*(invariant_mass - Higgs_mass)**2/(Higgs_gamma**2))
        return Higgs_gauss_factor    
    

    def Window(self, invariant_mass):
        mask = (invariant_mass >= self.window[0]) & (invariant_mass <= self.window[1])
        return mask

    #This is experimental part and by default not used by main code

    def ptLikelihood(self, pTTauTau: np.ndarray, type: np.ndarray):

        mask1 = (np.abs(pTTauTau)<0.5)

        if type == 0:
            pT1 = self.leg1P4[:, 0][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
            pT2 = self.leg2P4[:, 0][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
        elif type == 1:
            pT1 = self.leg1P4[:, 1][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
            pT2 = self.leg2P4[:, 1][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
        elif type == 2:
            pT1 = self.leg1P4[:, 2][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
            pT2 = self.leg2P4[:, 2][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))

        x1Min = np.minimum(1.0, self.mVisOverTauSquare1)[:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
        x2Min = np.minimum(1.0, self.mVisOverTauSquare2)[:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))

        x1Max = np.ones(pTTauTau.shape)
        x2Max = np.ones(pTTauTau.shape)

        a_x2 = x1Min *pT2/(x1Min*pTTauTau - pT1)
        b_x2 = x1Max*pT2/(x1Max*pTTauTau - pT1)

        x1_singularity = pT1/pTTauTau
        x2_vs_x1_singularity = (x1_singularity>0.0) & (x1_singularity<1.0)

        momentum_sign = (-pT2*pT1<0)

        x2Min = np.where(momentum_sign, np.maximum(x2Min, b_x2), x2Min)
        x2Max = np.where(momentum_sign, np.minimum(x2Max, a_x2), x2Max)
        x2Max = np.where((momentum_sign) & (x2_vs_x1_singularity) & (x2Max<0), 1.0, x2Max)
        x2Min = np.where(~momentum_sign, np.maximum(x2Min, a_x2), x2Min)
        x2Max = np.where(~momentum_sign, np.minimum(x2Max, b_x2), x2Max)
        x2Max = np.where((~momentum_sign) & (x2_vs_x1_singularity) & (x2Max<0), 1.0, x2Max)

        x2Min[x2Min<0] = 0.0
        
        mask2 = (x2Min > x2Max)

        HadDecay1 = np.broadcast_to((self.leg1DecayType != 1)[:, np.newaxis], pTTauTau.shape)
        HadDecay2 = np.broadcast_to((self.leg2DecayType != 1)[:, np.newaxis], pTTauTau.shape)
        
        mNuNuIntegral = np.zeros((pTTauTau.shape))
        x2 = np.minimum(1.0, x2Max)

        term1 = pT2 - pTTauTau*x2
        log_term1 = np.log(np.abs(term1))

        integralMax = pT1*(pTTauTau*x2 + pT2**2/term1 + 2*pT2*log_term1)/pTTauTau**3

        ###MOST CONSUMING PART 1###

        mNuNuIntegral += HadDecay1 * (-pT1**2*(2*pTTauTau*x2+pT2**2*(5*pT2-6*pTTauTau*x2)/term1**2 + 6*pT2*log_term1)/(2*pTTauTau**4))
        mNuNuIntegral += HadDecay2 * (-pT1/(2*pTTauTau**5)*(2*pT2*pTTauTau*(-3*pT1 + 2*pTTauTau)*x2 + pTTauTau**2*(-pT1 + pTTauTau)*x2**2 + (pT2**4*pT1)/term1**2 + 2*pT2**3*(-4*pT1 + pTTauTau)/term1 + 6*pT2**2*(-2*pT1 + pTTauTau)*log_term1))

        integralMax += mNuNuIntegral

        ###END OF MOST CONSUMING PART 1###

        mNuNuIntegral = np.zeros((pTTauTau.shape))

        x2 = x2Min
        term2 = pT2 - pTTauTau*x2
        log_term2 = np.log(np.abs(term2))

        integralMin = pT1*(pTTauTau*x2+pT2**2/term2+2*pT2*log_term2)/pTTauTau**3

        ###MOST CONSUMING PART 2###
        
        mNuNuIntegral += HadDecay1 * (-pT1**2*(2*pTTauTau*x2+pT2**2*(5*pT2-6*pTTauTau*x2)/term2**2+6*pT2*log_term2)/(2*pTTauTau**4))
        mNuNuIntegral += HadDecay2 * (-pT1/(2*pTTauTau**5)*(2*pT2*pTTauTau*(-3*pT1 + 2*pTTauTau)*x2 + pTTauTau**2*(-pT1 + pTTauTau)*x2**2 + (pT2**4*pT1)/term2**2 + 2*pT2**3*(-4*pT1 + pTTauTau)/term2 + 6*pT2**2*(-2*pT1 + pTTauTau)*log_term2))
        
        integralMin += mNuNuIntegral

        ###END OF MOST CONSUMING PART 2###

        value = integralMax - integralMin

        value[mask1 | mask2] = 0.0

        #value*=1E4

        return np.abs(value)
    
    def metTF(self, metP4: np.ndarray, nuP4: np.ndarray, covMET: np.ndarray) -> np.ndarray:
        aMETx = metP4[..., 0]
        aMETy = metP4[..., 1]

        covDET = np.linalg.det(covMET)
        mask = covDET < 1E-10
        covDET[mask] = 1.0

        constMET = 1/2/np.pi/np.sqrt(covDET)
        residualX = aMETx[:, np.newaxis] - nuP4[:, :, 0]
        residualY = aMETy[:, np.newaxis] - nuP4[:, :, 1]

        #covMET 1 coordinate responds to X and 0 coordinate to Y
        pull2 = residualX*(covMET[:, np.newaxis, 1, 1]*residualX - covMET[:, np.newaxis, 0, 1]*residualY) + residualY*(-covMET[:, np.newaxis, 1, 0]*residualX + covMET[:, np.newaxis, 0, 0]*residualY)
        pull2 /= covDET[:, np.newaxis]
        
        pull2[np.broadcast_to(mask[:, np.newaxis], pull2.shape)] = 0.0
        return constMET[:, np.newaxis]*np.exp(-0.5*pull2)
    
    def value(self, x: np.ndarray):
        
        x1Min = np.minimum(1.0, self.mVisOverTauSquare1)
        x2Min = np.minimum(1.0, self.mVisOverTauSquare2)

        mask = (x[:, 0] < x1Min[:, np.newaxis]) | (x[:, 1] < x2Min[:, np.newaxis])
        
        testP4 = self.leg1P4[:, np.newaxis, :] / x[:, 0][:, np.newaxis] + self.leg2P4[:, np.newaxis, :] / x[:, 1][:, np.newaxis]

        testMET = testP4 - self.leg1P4[:, np.newaxis, :] - self.leg2P4[:, np.newaxis, :]

        value = np.full(testMET.shape[:2], -1.0) #Negative likelihood
        
        if self.enable_MET:
            value *= self.metTF(self.recoMET, testMET, self.covMET)
        if self.enable_mass:
            value *= self.massLikelihood(InvariantMass(testP4))

        #Experimental components
        #Not  introduced yet in official version
        if self.enable_px:
            value *= self.ptLikelihood(testP4[:, :, 0], 0)
        if self.enable_py:
            value *= self.ptLikelihood(testP4[:, :, 1], 1)
        if self.enable_BW:
            value *= self.BreitWigner(InvariantMass(testP4))
        if self.enable_gauss:
            value *= self.Gauss(InvariantMass(testP4))
        if self.enable_window:
            value *= self.Window(InvariantMass(testP4))
            
        value[mask] = 0.000001

        return value

class FastMTT(Likelihood):
    def __init__(self,
                 enable_BW = False,
                 enable_window = False,
                 calculate_uncertainties = False):
        self.myLikelihood = Likelihood(enable_BW = False,
                                       enable_window = False)
        self.BestLikelihood = 0.0
        self.BestX = np.array([0.0, 0.0])
        self.bestP4 = 0.0
        self.tau1P4 = 0.0
        self.tau2P4 = 0.0
        self.mass = 0.0

        #New component to calculate uncertainty
        #It produces long tails, but apart from that calculates uncertainties event by event quite ok ~ after some cuts results are aprox. Gaussian
        #A bit time consuming -- doubles the time of calculation -- so it is disabled by default

        self.CalculateUncertainties = calculate_uncertainties
        self.one_sigma = 0.0

        #Number of event for which likelihood plot will be shown.
        #-1 = no plot
        self.WhichLikelihoodPlot = -1

        return
    
    def run(self, measuredTauLeptons, measuredMETx, measuredMETy, covMET) -> np.ndarray:

        start_real_time = time.time()
        start_cpu_time = time.process_time()

        ##############################################
                            #RUN
        ##############################################

        if np.shape(measuredTauLeptons)[1] != 2:
            print(f"Number of MeasuredTauLepton is {len(measuredTauLeptons)}. A user shouls pass exactly two leptons.\n")
            return

        metLenght = np.sqrt(measuredMETx**2 + measuredMETy**2)
        aMET = np.array([measuredMETx, measuredMETy, np.zeros(np.shape(measuredMETx)), metLenght]).T

        aLepton1 = measuredTauLeptons[:, 0]
        aLepton2 = measuredTauLeptons[:, 1]

        self.p4_Lepton1 = self.get_p4(aLepton1)
        self.p4_Lepton2 = self.get_p4(aLepton2)

        self.Lepton1 = self.p4_Lepton1
        self.Lepton2 = self.p4_Lepton2

        aLepton1 = self.modify_lepton_mass(aLepton1)
        aLepton2 = self.modify_lepton_mass(aLepton2)

        self.myLikelihood.mvisleg1 = aLepton1[:, 4]
        self.myLikelihood.mvisleg2 = aLepton2[:, 4]

        #setMETinputs
        self.myLikelihood.recoMET = aMET
        self.myLikelihood.covMET = covMET

        self.myLikelihood.setLeptonInputs(self.p4_Lepton1, self.p4_Lepton2, aLepton1[:, 0], aLepton2[:, 0], aLepton1[:, 5], aLepton2[:, 5])

        self.scan()
        
        self.tau1P4 = self.p4_Lepton1*(1/self.BestX[:, np.newaxis, 0])
        self.tau2P4 = self.p4_Lepton2*(1/self.BestX[:, np.newaxis, 1])
        self.bestP4 = self.tau1P4 + self.tau2P4
        self.mass = InvariantMass(self.bestP4)

        #if self.myLikelihood.enable_window:
        #    self.mass = self.mass[(self.mass >= self.myLikelihood.window[0]) & (self.mass <= self.myLikelihood.window[1])]

        ##if self.CalculateUncertainties == True:
        ##    self.propagate_uncertainties()

        ##############################################

        #Time calculation part:
        end_real_time = time.time()
        end_cpu_time = time.process_time()
        
        real_time_elapsed = round(end_real_time - start_real_time, 3)
        cpu_time_elapsed = round(end_cpu_time - start_cpu_time, 3)

        logger.info(f"Time elapsed : CPU - {cpu_time_elapsed} & Real - {real_time_elapsed} seconds")
    
    #lepton[0]: decay_type:
    #1 - TauToHad
    #2 - TauToElec
    #3 - TauToMu

    #lepton[1]: pt
    #lepton[2]: eta
    #lepton[3]: phi
    #lepton[4]: mass
    #lepton[5]: hadron decay mode

    def get_p4(self, lepton: np.ndarray):
        p = lepton[:, 1] * np.cosh(lepton[:, 2])
        px = lepton[:, 1] * np.cos(lepton[:, 3])
        py = lepton[:, 1] * np.sin(lepton[:, 3])
        pz = lepton[:, 1] * np.sinh(lepton[:, 2])
        energy = np.sqrt(p**2 + lepton[:, 4])
        return np.array([px, py, pz, energy]).T
    
    def modify_lepton_mass(self, aLepton1, electron_mass=ElectronMass, muon_mass=MuonMass, pion_mass=ChargedPionMass):

        # 1) Set electron inv mass to electron
        aLepton1[aLepton1[:, 0] == 2, 4] = electron_mass

        # 2) Set muon inv mass to muon
        aLepton1[aLepton1[:, 0] == 3, 4] = muon_mass

        # 3) If decay is hadronic, we set minimal and maximal mass, depending on decay mode
        mask_type1 = aLepton1[:, 0] == 1

        # a) If hadronic decay mode == -1, then min mass = pion_mass and max mass = 1.5
        mask_a = (mask_type1 & (aLepton1[:, 5] == -1) & (aLepton1[:, 4] < pion_mass))
        aLepton1[mask_a, 4] = pion_mass

        mask_b = (mask_type1 & (aLepton1[:, 5] == -1) & (aLepton1[:, 4] > 1.5))
        aLepton1[mask_b, 4] = 1.5

        # c) If hadronic decay mode == 0, then mass = pion mass
        mask_c = (mask_type1 & (aLepton1[:, 5] == 0))
        aLepton1[mask_c, 4] = pion_mass

        # d) If hadronic decay mode has another value, we set min mass to 0.3 and max mass to 1.5
        mask_d = (mask_type1 & (aLepton1[:, 5] != -1) & (aLepton1[:, 5] != 0) & (aLepton1[:, 4] < 0.3))
        aLepton1[mask_d, 4] = 0.3

        mask_e = (mask_type1 & (aLepton1[:, 5] != -1) & (aLepton1[:, 5] != 0) & (aLepton1[:, 4] > 1.5))
        aLepton1[mask_e, 4] = 1.5

        return aLepton1
    
    def scan(self):
        
        nGridPoints = 100
        gridFactor = 1.0/nGridPoints

        X1 = np.arange(1, nGridPoints) * gridFactor
        X2 = np.arange(1, nGridPoints) * gridFactor

        # Cartesian product
        self.pairs = np.column_stack((np.repeat(X1, len(X2)),
                                 np.tile(X2, len(X1))))
        
        self.lh = self.myLikelihood.value(self.pairs)

        minimum = np.argmin(self.lh, axis=1)

        self.BestX = self.pairs[minimum]
        self.BestLikelihood = self.lh[np.arange(self.lh.shape[0]), minimum]

        ### USER INTERFACE AND ADDITIONAL COMPONENTS ###

        chi_square = 2.3

        #Plotting likelihoods
        #if self.WhichLikelihoodPlot != -1:
        #    self.plot_likelihood(X1, X2, event_number = self.WhichLikelihoodPlot, threshold=self.BestLikelihood[self.WhichLikelihoodPlot]/np.exp(chi_square/2))

        if self.CalculateUncertainties == True:
            self.contour_uncertainties(X1, X2, chi_square)

        #Code for minimalizing function with scipy:

        '''initial_guess = np.array([0.5, 0.5])
        result = minimize(self.myLikelihood.value, initial_guess, method='BFGS')
        self.BestX = result.x
        self.BestLikelihood = result.fun'''
        
        #Faster than grid search in pure python
        #Slower than grid search in numpy with vectorization and broadcasting
        #Potentially one can replace it with jax and/or numba

        return

    """
    def plot_likelihood(self, X1, X2, event_number=0, threshold=None):

        print("Threshold: ", threshold)
        nGridPoints = np.shape(X1)[0]

        lh_grid = self.lh[event_number, :].reshape(nGridPoints, nGridPoints)

        plt.figure(figsize=(8, 6))

        # Main heatmap
        img = plt.imshow(-lh_grid, origin='lower', extent=(X1[0], X1[-1], X2[0], X2[-1]), 
                        cmap='viridis_r', interpolation='nearest')

        # Unphysical region excluded
        unphysical_region_mask = (lh_grid == 0.000001)
        plt.contourf(X1, X2, unphysical_region_mask, levels=[0.5, 1], colors='red', alpha=0.6)
        contour_proxy = plt.Line2D([0], [0], linestyle="none", marker="s", markersize=10, markerfacecolor="red", alpha=0.6)
        plt.legend([contour_proxy], [r"Kinematic limit $(\frac{m_i}{m_\tau})^2<x_i$"], loc="upper left")

        # Maximum likelihood point
        maximum_likelihood = self.BestX[event_number]
        plt.scatter(maximum_likelihood[1], maximum_likelihood[0], color='red', marker='x', s=100, label="Maximum likelihood", linewidth=2.5)

        # Contour for uncertainties (threshold = 2.3 for 1 sigma)
        if threshold is not None:
            plt.contour(X1, X2, lh_grid, levels=[threshold], colors='blue', linewidths=2, linestyles='--')

        # Finalizing
        plt.colorbar(img, label='Likelihood')

        plt.xlabel(r"$X_1$")
        plt.ylabel(r"$X_2$")
        plt.title("2D Heatmap of Likelihood Function")

        params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (10, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
        plt.rcParams.update(params)

        file_path = f"images/fastMTT/likelihood_{self.WhichLikelihoodPlot}_event.png"
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        plt.savefig(file_path, dpi=300)
        plt.close()
    """
    def evaluate_mass(self, x):
        tau1P4 = self.p4_Lepton1[:, np.newaxis, np.newaxis, :]*(1/x[np.newaxis, :, :, np.newaxis, 0])
        tau2P4 = self.p4_Lepton2[:, np.newaxis, np.newaxis, :]*(1/x[np.newaxis, :, :, np.newaxis, 1])
        bestP4 = tau1P4 + tau2P4
        mass = InvariantMass(bestP4)
        return mass
    
    def contour_uncertainties(self, X1, X2, chi_square = 2.3):
        threshold = self.BestLikelihood/np.exp(chi_square/2)

        nGridPoints = np.shape(X1)[0]
        nEvents = np.shape(self.lh)[0]

        lh_grid = self.lh.reshape(nEvents, nGridPoints, nGridPoints)
        pairs = self.pairs.reshape(nGridPoints, nGridPoints, 2)
        mask = (lh_grid < threshold[:, np.newaxis, np.newaxis])

        up = np.roll(mask, shift=-1, axis=1)
        down = np.roll(mask, shift=1, axis=1)
        left = np.roll(mask, shift=-1, axis=2)
        right = np.roll(mask, shift=1, axis=2)
        boundary_mask = mask & ~(up & down & left & right)
        
        masses = self.evaluate_mass(pairs)
        masses = np.where(boundary_mask, masses, np.nan)
        max_masses = np.nanmax(masses, axis=(1, 2))
        min_masses = np.nanmin(masses, axis=(1, 2))

        arbitrary_numerical_factor = 2.5

        self.one_sigma = (max_masses - min_masses)/2*arbitrary_numerical_factor
