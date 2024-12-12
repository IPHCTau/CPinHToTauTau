import os
from columnflow.util import maybe_import
from typing import Optional

np     = maybe_import("numpy")
ak     = maybe_import("awkward")
coffea = maybe_import("coffea")

from httcp.production.TComplex import TComplex
from columnflow.columnar_util import EMPTY_FLOAT, Route, flat_np_view, layout_ak_array




class PolarimetricA1:
    def __init__(
            self,
            p4_tau: ak.Array, 
            p4_os_pi: ak.Array, p4_ss1_pi: ak.Array, p4_ss2_pi: ak.Array,
            taucharge: ak.Array,
            decayChannel: Optional[str]="k3ChargedPi"
    ) -> None:
        """
          Calculate Polarimetric vector for tau to a1 decay
          args:
            p4_tau, p4_os_pi, p4_ss1_pi, p4_ss2_pi : awkward array of lorentz vectors
            taucharge
            decayChannel
          returns:
            None
        """
        print("Configure PolarimetricA1 : Start")

        self.P   = p4_tau
        self.p1  = p4_os_pi
        self.p2  = p4_ss1_pi
        self.p3  = p4_ss2_pi
        self.taucharge = taucharge
        self.decayChannel = decayChannel
        
        self.numResonances  = 7
        # define mass of tau lepton, charged and neutral pion;
        # values are taken from Prog. Theor. Exp. Phys. 2022 (2022) 083C01 (PDG)
        self.m_tau          = 1.7769   # [GeV]
        self.m_chargedPi    = 0.139570 # [GeV]
        self.m_neutralPi    = 0.134977 # [GeV]
        
        # define mass and width of a1(1260) meson;
        # values are taken from column "nominal fit" in Table VI of Phys.Rev.D 61 (2000) 012002
        self.m0_a1          = 1.331    # [GeV]
        self.Gamma0_a1      = 0.814    # [GeV]

        # define parameters specifying "running" of a1 width;
        # the values of Gamma_a1 as function of s have been taken from Fig. 9 (b) of Phys.Rev.D 61 (2000) 012002
        Gamma_a1_vs_s = [
            ( 0.00, 0.000 ), ( 0.20, 0.000 ), ( 0.40, 0.000 ), ( 0.50, 0.005 ), ( 0.60, 0.020 ),
            ( 0.65, 0.040 ), ( 0.70, 0.055 ), ( 0.75, 0.075 ), ( 0.80, 0.110 ), ( 0.85, 0.160 ),
            ( 0.90, 0.205 ), ( 0.95, 0.250 ), ( 1.00, 0.295 ), ( 1.05, 0.340 ), ( 1.10, 0.375 ),
            ( 1.15, 0.410 ), ( 1.20, 0.450 ), ( 1.25, 0.475 ), ( 1.30, 0.515 ), ( 1.35, 0.555 ),
            ( 1.40, 0.595 ), ( 1.45, 0.630 ), ( 1.50, 0.660 ), ( 1.55, 0.690 ), ( 1.60, 0.720 ),
            ( 1.65, 0.750 ), ( 1.70, 0.780 ), ( 1.75, 0.815 ), ( 1.80, 0.845 ), ( 1.85, 0.875 ),
            ( 1.90, 0.905 ), ( 1.93, 0.930 ), ( 1.95, 0.980 ), ( 2.00, 1.060 ), ( 2.05, 1.125 ),
            ( 2.10, 1.185 ), ( 2.15, 1.245 ), ( 2.20, 1.300 ), ( 2.25, 1.355 ), ( 2.30, 1.415 ),
            ( 2.35, 1.470 ), ( 2.37, 1.485 ), ( 2.40, 1.520 ), ( 2.45, 1.575 ), ( 2.50, 1.640 ),
            ( 2.55, 1.705 ), ( 2.60, 1.765 ), ( 2.65, 1.835 ), ( 2.70, 1.900 ), ( 2.75, 1.970 ),
            ( 2.80, 2.050 ), ( 2.85, 2.130 ), ( 2.90, 2.205 ), ( 2.95, 2.285 ), ( 3.00, 2.380 ),
            ( 3.05, 2.470 ), ( 3.10, 2.570 ), ( 3.15, 2.690 ) 
        ]
        self.Gamma_a1_vs_s = [list(tup) for tup in Gamma_a1_vs_s]

        # define masses and widths of intermediate rho(770), rho(1450), f2(1270), sigma, f0(1370) resonances;
        # values are taken from Table I of Phys.Rev.D 61 (2000) 012002
        self.m0_rho770      = 0.774    # [GeV]
        self.Gamma0_rho770  = 0.149    # [GeV]
        self.m0_rho1450     = 1.370    # [GeV]
        self.Gamma0_rho1450 = 0.386    # [GeV]
        self.m0_f2          = 1.275    # [GeV]
        self.Gamma0_f2      = 0.185    # [GeV]
        self.m0_sigma       = 0.860    # [GeV]
        self.Gamma0_sigma   = 0.880    # [GeV]
        self.m0_f0          = 1.186    # [GeV]
        self.Gamma0_f0      = 0.350    # [GeV]

        # define coefficients specifying the contribution of meson resonances to the hadronic current J
        # values are taken from Table III of Phys.Rev.D 61 (2000) 012002
        self.beta_moduli = [  1.00,  0.12,  0.37,  0.87,  0.71,  2.10,  0.77 ]
        self.beta_phases = [  0.00,  0.99, -0.15,  0.53,  0.56,  0.23, -0.54 ]
        assert len(self.beta_moduli) == self.numResonances and len(self.beta_phases) == self.numResonances

        self.beta = [TComplex(self.beta_moduli[i], self.beta_phases[i]*np.pi, True) for i in range(self.numResonances)]

        # metric tensor g^{mu,nu}
        self.g = np.zeros((4,4), dtype=TComplex)
        for mu in range(4):
            for nu in range(4):
                self.g[mu][nu] = self.get_g(mu, nu)
        #self.g = np.array([[-1.0, -1.0, -1.0, -1.0],
        #                   [ 0.0, -1.0,  0.0, -1.0],
        #                   [ 0.0,  0.0, -1.0, -1.0],
        #                   [ 0.0,  0.0,  0.0,  1.0]])
        
        
    def get_g(self, mu: int, nu: int) -> int:
        if (mu == 0 and nu == 0) or (mu == 1 and nu == 1) or (mu == 2 and nu == 2): return -1.0
        elif mu == 3 and nu == 3:                                                   return 1.0
        else:                                                                       return 0.0

        
    def Boost(self,p,frame):
        boostvec = frame.boostvec
        return p.boost(boostvec.negative())
        

    def Setup(self):
        pass

    
    def PVC(self) -> ak.Array:
        #dummy = ak.zeros_like(self.P.pt)
        #P = setp4("PtEtaPhiMLorentzVector", dummy, dummy, dummy, self.P.mass)
        #P  = self.Boost(self.P, self.P)
        P  = self.P
        p1 = self.p1
        p2 = self.p2
        p3 = self.p3

        # saving original dimension
        numP = ak.num(P, axis=1)
        mask = numP > 0
        P = P[mask]
        p1 = p1[mask]
        p2 = p2[mask]
        p3 = p3[mask]
        tauc = self.taucharge[mask]
        
        #_layout = ak.num(P,  axis=1)[:,None]
        
        #from IPython import embed; embed()
        
        #p1 = self.Boost(self.p1, self.P)
        #p2 = self.Boost(self.p2, self.P)
        #p3 = self.Boost(self.p3, self.P)
        
        N = P - (p1 + p2 + p3)

        J = self.comp_J(p1, p2, p3)
        #print(f" >>> J: {J[0]}, {J[1]}, {J[2]}, {J[3]}")
        #print(f"J[0]: {J[0].Re()}, {J[0].Im()}, {type(J[0].Im())}")
        #print(f"J[1]: {J[1].Re()}, {J[1].Im()}, {type(J[1].Im())}")
        
        Pi  = self.comp_Pi(J, N)
        print(f" >>> Pi: {Pi}")
        
        # charge
        Pi5 = self.comp_Pi5(J, N)
        print(f" >>> Pi5: {Pi5}")
        
        # CV: Standard Model value, cf. text following Eq. (3.15) in Comput.Phys.Commun. 64 (1991) 275
        gammaVA = 1.0

        # CV: sign of terms proportional to gammaVA differs for tau+ and tau-,
        #     cf. text following Eq. (3.16) in Comput.Phys.Commun. 64 (1991) 275
        sign = np.where(tauc > 0, -1, 1)
        #if    self.taucharge == +1: sign = -1.
        #elif  self.taucharge == -1: sign = +1.
        #else: assert False

        print(f"sign: {sign}")
        
        omega = P.dot(Pi - sign*gammaVA*Pi5)
        print(f"omega: {omega}, {type(omega)}")

        print(f"1./(omega*self.m_tau): {1./(omega*self.m_tau)}")
        print(f"np.power(self.m_tau, 2)*(Pi5 - sign*gammaVA*Pi): {np.power(self.m_tau, 2)*(Pi5 - sign*gammaVA*Pi)}")
        print(f"P.dot(Pi5 - sign*gammaVA*Pi)*P: {P.dot(Pi5 - sign*gammaVA*Pi)*P}")
        
        #P = setp4("LorentzVector", P.px, P.py, P.pz, P.energy)
        P = ak.zip(
            {
                "x": P.px,
                "y": P.py,
                "z": P.pz,
                "t": P.energy,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )

        
        #H = (1./(omega*P.mass))*(np.power(P.mass, 2)*(Pi5 - sign*gammaVA*Pi) - P.dot(Pi5 - sign*gammaVA*Pi)*P)
        H = (1./(omega*self.m_tau))*(np.power(self.m_tau, 2)*(Pi5 - sign*gammaVA*Pi) - P.dot(Pi5 - sign*gammaVA*Pi)*P)
        print(f"H: {H}, {type(H)}")
        
        #retVal = H.pvec.unit
        #retVal = H.pvec

        #return retVal

        #H = layout_ak_array(np.flatten(H), _layout)


        fH = ak.flatten(H)
        masked_numP = numP.mask[numP == 1]
        flat_index = np.cumsum(masked_numP) - 1
        flat_result = fH[flat_index]
        H = flat_result[:,None]
        H = ak.drop_none(H)
        
        return H

        
    def comp_J(self, p1: ak.Array, p2: ak.Array, p3: ak.Array) -> np.array(4, dtype=TComplex):
        """
        Args:
          p1, p2, p3: awkward array of Lorentz Vectors

        Return:
          Awkward array of TComplex objects
        """
        print(" --- comp_J --- ")
        print(f" p1: {p1}, {p1.type}")
        q1 = p2 - p3
        q2 = p3 - p1
        q3 = p1 - p2

        h1 = p2 + p3
        Q1 = h1 - p1
        s1 = h1.mass2
        h2 = p1 + p3
        Q2 = h2 - p2
        s2 = h2.mass2
        h3 = p1 + p2
        Q3 = h3 - p3
        s3 = h3.mass2
        
        m1, m2, m3 = 0., 0., 0.
        if self.decayChannel == "k3ChargedPi":
            m1 = self.m_chargedPi
            m2 = self.m_chargedPi
            m3 = self.m_chargedPi
        elif self.decayChannel == "kChargedPi2NeutralPi":
            m1 = self.m_neutralPi
            m2 = self.m_neutralPi
            m3 = self.m_chargedPi
        else:
            raise RuntimeError(f"Error in <PolarimetricVectorTau2a1::comp_J>: Invalid parameter 'decayChannel' = {self.decayChannel}")

  
        a = p1 + p2 + p3
        s = a.mass2

        #from IPython import embed; embed()
        
        T = np.zeros((4,4), dtype=TComplex)
        for mu in range(4):
            for nu in range(4):
                # CV: all vectors without explicit mu indices are assumed to have the mu index as subscript,
                #     hence multiplication with the product of metric tensors g^{mu,mu}*g^{nu,nu} is required to transform
                #     the return value of the functions get_component(a, mu)*get_component(a, nu) into the expression a^{mu}*a^{nu}
                #     when computing T^{mu,nu} = g^{mu,nu} - a^{mu}a^{nu}/a^2,
                #     as described in text following Eq. (A2) in Phys.Rev.D 61 (2000) 012002
                comp1np = ak.to_numpy(self.get_component(a, mu))
                comp2np = ak.to_numpy(self.get_component(a, nu))
                snp     = ak.to_numpy(s)
                #T[mu][nu] = TComplex(self.g[mu][nu] - self.g[mu][mu]*self.g[nu][nu]*self.get_component(a, mu)*self.get_component(a, nu)/s)
                #T[mu][nu] = TComplex(self.g[mu][nu] - self.g[mu][mu]*self.g[nu][nu]*self.get_component(a, mu)*self.get_component(a, nu)/ak.to_numpy(s))
                T[mu][nu] = TComplex(self.g[mu][nu] - self.g[mu][mu]*self.g[nu][nu]*comp1np*comp2np/snp)
                

        #print(f"T: {T}")
        #print(type(T[0][2].Re()))
        #print(type(T[0][2].Im()))
        #for mu in range(4):
        #    for nu in range(4):
        #        print(f"T_Re: {T[mu][nu].Re()}")
        #        print(f"T_Im: {T[mu][nu].Im()}")
        # compute amplitudes for individual resonances according to Eq. (A3) in Phys.Rev.D 61 (2000) 012002
        # Note: all the factors F_R in Eq. (A3) are equal to one, as the nominal fit assumes the a1 size parameter R to be zero
        #j = [np.zeros((4,snp.shape[0],1), dtype=TComplex) for _ in range(self.numResonances)]


        j = [np.zeros(4, dtype=TComplex) for _ in range(self.numResonances)]
        cq1 = self.convert_to_cLorentzVector(q1)
        cq2 = self.convert_to_cLorentzVector(q2)

        # ----------------------------------------------------------------------------------------------------------------- #
        # ----------------------------------------------------- j [0] ----------------------------------------------------- #
        # ----------------------------------------------------------------------------------------------------------------- #
        print(" ----------------------------------------------------- j [0] ----------------------------------------------------- ")
        #j[0] = T@(self.BreitWigner(self.m0_rho770,  self.Gamma0_rho770,  s1, m2, m3, 1)*cq1 - self.BreitWigner(self.m0_rho770,  self.Gamma0_rho770,  s2, m1, m3, 1)*cq2)
        #j[1] = T@(self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s1, m2, m3, 1)*cq1 - self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s2, m1, m3, 1)*cq2)
        j01   = self.get_compJ_comps(self.BreitWigner(self.m0_rho770,  self.Gamma0_rho770,  s1, m2, m3, 1), cq1)
        j02   = self.get_compJ_comps(self.BreitWigner(self.m0_rho770,  self.Gamma0_rho770,  s2, m1, m3, 1), cq2)
        j0    = j01 - j02
        print(f"j0: {j0}")
        j[0]  = T@j0
        #j[0] = self.get_mult_comps(T, j0)
        print(f"j[0]: {j[0]}, {j[0].shape}")

        # ----------------------------------------------------------------------------------------------------------------- #
        # ----------------------------------------------------- j [1] ----------------------------------------------------- #
        # ----------------------------------------------------------------------------------------------------------------- #
        print(" ----------------------------------------------------- j [1] ----------------------------------------------------- ")
        j11   = self.get_compJ_comps(self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s1, m2, m3, 1), cq1)
        j12   = self.get_compJ_comps(self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s2, m1, m3, 1), cq2)
        j1    = j11 - j12
        #j1   = self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s1, m2, m3, 1)*cq1 - self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s2, m1, m3, 1)*cq2
        #j[1] = np.sum(T * j1.T, axis=1)
        j[1] = T@j1
        print(f"j[1]: {j[1]}, {j[1].shape}")
        
        # ----------------------------------------------------------------------------------------------------------------- #
        # ----------------------------------------------------- j [2] ----------------------------------------------------- #
        # ----------------------------------------------------------------------------------------------------------------- #
        print(" ----------------------------------------------------- j [2] ----------------------------------------------------- ")
        aXq1 = ak.to_numpy(a.dot(q1))
        print(f"aXq1: {aXq1}")
        cQ1 = self.convert_to_cLorentzVector(Q1)
        
        aXq2 = ak.to_numpy(a.dot(q2))
        print(f"aXq2: {aXq2}")        
        cQ2 = self.convert_to_cLorentzVector(Q2)
        
        #j[2] = T@(aXq1*self.BreitWigner(self.m0_rho770,  self.Gamma0_rho770,  s1, m2, m3, 1)*cQ1 - aXq2*self.BreitWigner(self.m0_rho770,  self.Gamma0_rho770,  s2, m1, m3, 1)*cQ2)
        #j[3] = T@(aXq1*self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s1, m2, m3, 1)*cQ1 - aXq2*self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s2, m1, m3, 1)*cQ2)
        j21   = self.get_compJ_comps(self.BreitWigner(self.m0_rho770, self.Gamma0_rho770,  s1, m2, m3, 1), cQ1)
        j21   = self.get_mult_arrs(j21, aXq1)

        j22   = self.get_compJ_comps(self.BreitWigner(self.m0_rho770,  self.Gamma0_rho770,  s2, m1, m3, 1), cQ2) 
        j22   = self.get_mult_arrs(j22, aXq2)

        j2    = j21 - j22
        j[2]  = T@j2
        print(f"j[2]: {j[2]}, {j[2].shape}")

        # ----------------------------------------------------------------------------------------------------------------- #
        # ----------------------------------------------------- j [3] ----------------------------------------------------- #
        # ----------------------------------------------------------------------------------------------------------------- #
        print(" ----------------------------------------------------- j [3] ----------------------------------------------------- ")
        j31  = self.get_compJ_comps(self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s1, m2, m3, 1), cQ1)
        j31  = self.get_mult_arrs(j31, aXq1)

        j32  = self.get_compJ_comps(self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s2, m1, m3, 1), cQ2)
        j32  = self.get_mult_arrs(j32, aXq2)
        
        #j3   = aXq1*self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s1, m2, m3, 1)*cQ1 - aXq2*self.BreitWigner(self.m0_rho1450, self.Gamma0_rho1450, s2, m1, m3, 1)*cQ2
        j3   = j31 - j32
        #j[3] = np.sum(T * j3.T, axis=1)
        j[3] = T@j3
        print(f"j[3]: {j[3]}, {j[3].shape}")
        

        # ----------------------------------------------------------------------------------------------------------------- #
        # ----------------------------------------------------- j [4] ----------------------------------------------------- #
        # ----------------------------------------------------------------------------------------------------------------- #
        print(" ----------------------------------------------------- j [4] ----------------------------------------------------- ")
        aXq3 = ak.to_numpy(a.dot(q3))
        cq3 = self.convert_to_cLorentzVector(q3)
        q3Xq3 = ak.to_numpy(q3.mass2)
        ca = self.convert_to_cLorentzVector(a)
        h3Xa = ak.to_numpy(h3.dot(a))
        ch3 = self.convert_to_cLorentzVector(h3)

        #j[4] = T@(self.BreitWigner(self.m0_f2, self.Gamma0_f2, s3, m1, m2, 2)*(aXq3*cq3 - (q3Xq3/3.)*(ca - (h3Xa/s3)*ch3)))
        j4a  = self.BreitWigner(self.m0_f2, self.Gamma0_f2, s3, m1, m2, 2)
        j4b  = self.get_mult_arrs(cq3, aXq3) - self.get_mult_arrs((ca - self.get_mult_arrs(ch3, (h3Xa/s3))), (q3Xq3/3.))
        #j4   = self.BreitWigner(self.m0_f2, self.Gamma0_f2, s3, m1, m2, 2)*(aXq3*cq3 - (q3Xq3/3.)*(ca - (h3Xa/s3)*ch3))
        j4   = self.get_compJ_comps(j4a, j4b)
        #j[4] = np.sum(T * j4.T, axis=1)
        j[4] = T@j4
        print(f"j[4]: {j[4]}, {j[4].shape}")


        # ----------------------------------------------------------------------------------------------------------------- #
        # ----------------------------------------------------- j [5] ----------------------------------------------------- #
        # ----------------------------------------------------------------------------------------------------------------- #
        print(" ----------------------------------------------------- j [5] ----------------------------------------------------- ")
        cQ3  = self.convert_to_cLorentzVector(Q3)
        j5   = self.get_compJ_comps(self.BreitWigner(self.m0_sigma, self.Gamma0_sigma, s3, m1, m2, 0), cQ3)
        #j[5] = np.sum(T * j5.T, axis=1)[:,:,np.newaxis]
        j[5] = T@j5 
        print(f"j[5]: {j[5]}, {j[5].shape}")

        
        # ----------------------------------------------------------------------------------------------------------------- #
        # ----------------------------------------------------- j [6] ----------------------------------------------------- #
        # ----------------------------------------------------------------------------------------------------------------- #
        print(" ----------------------------------------------------- j [6] ----------------------------------------------------- ")
        j6   = self.get_compJ_comps(self.BreitWigner(self.m0_f0, self.Gamma0_f0, s3, m1, m2, 0), cQ3)
        #j[6] = np.sum(T * j6.T, axis=1)[:,:,np.newaxis]
        j[6] = T@j6
        print(f"j[6]: {j[6]}, {j[6].shape}")


        #print("sdcdsc", j[6] + j[5])

        #retVal = np.zeros((4, snp.shape[0], 1) dtype=TComplex)
        #retVal = np.zeros(4, dtype=TComplex)
        retVal = self.get_mult_arrs(j[0], self.beta[0])
        for idx in range(self.numResonances):
            #retVal = retVal + j[idx]*self.beta[idx]
            if idx == 0: continue
            retVal += self.get_mult_arrs(j[idx], self.beta[idx])
            
        #retVal *= self.BreitWigner_a1(s)
        retVal = self.get_mult_arrs(retVal, self.BreitWigner_a1(s))

        #for i in range(4):
        #    print(f"alscascnjascn {i}:\n{retVal[i].Re()}\n{retVal[i].Im()}\n")
        
        
        # CV: multiply with metric tensor in order to transform J^{mu} into J_{mu},
        #     as required to insert J^{mu} given by Eq. (A2) of Phys.Rev.D 61 (2000) 012002 into
        #     Eq. (3.15) of Comput.Phys.Commun. 64 (1991) 275
        #retVal = ((retVal.T)@(self.g.T)).T
        retVal = self.g@retVal
        #for i in range(4):
        #    print(f"{i}:\n{retVal[i].Re()}\n{retVal[i].Im()}\n")
        return retVal



    
    def get_mult_arrs(self, TCarr, arr):
        """
        TCarr :  np.array(4, dtype=TComplex) -> each entry: [[], [], [], ..., []]
        arr: [[], [], [], ..., []]
        """
        retval = np.zeros(4, dtype=TComplex)
        for i in range(TCarr.shape[0]):
            retval[i] = TCarr[i]*arr
        return retval

    
    def get_sum_arrs(self, TCarr, arr):
        """
        TCarr :  np.array(4, dtype=TComplex) -> each entry: [[], [], [], ..., []]
        arr: [[], [], [], ..., []]
        """
        retval = np.zeros(4, dtype=TComplex)
        for i in range(TCarr.shape[0]):
            retval[i] = TCarr[i] + arr
        return retval

    
    
    def get_compJ_comps(self, BW, cLV):
        """
        cLV: numpy array of TComplex
             [TCx, TCy, TCz, TCt] (4, ) : TCx -> [[21.1], [22.3], [65.1], ..., [12.5]] (nEvents, 1)
        BW : numpy array of TComplex
             [[TC], [TC], [TC], [TC], ...., [TC]]
        """
        print(" --- get_compJ_comps --- ")
        comp = np.zeros(4, dtype=TComplex)
        for i in range(cLV.shape[0]):
            comp[i] = BW * cLV[i]
        print(comp)
        return comp

    
    def get_mult_comps(self, T, cLV):
        """
        cLV: numpy array of TComplex
             [TCx, TCy, TCz, TCt] (4, ) : TCx -> [[21.1], [22.3], [65.1], ..., [12.5]] (nEvents, 1)
        T  : numpy array of TComplex (4,4)
             each entry: TC -> [[10.1], [2.5], [3.4], ...., [12.7]]
        """
        import copy
        print(" --- get_compJ_comps --- ")

        TcLV = np.zeros(4, dtype=TComplex)
        for i in range(4):
            print(i)
            TcLVi = np.zeros(4, dtype=TComplex)
            for j in range(4):
                print(f"j: {j}")
                print(T[i][j])
                print(T[i][j].Re(), T[i][j].Im())
                print(cLV[i])
                print(cLV[i].Re(), cLV[i].Im())
                TcLV0[i] = T[i][j]*cLV[i]
                print(TcLV0[i])
            TcLVi = np.sum(TcLV, axis=1)
            print(TcLVi)
            TcLV[i] = copy.deepcopy(TcLVi)

        return TcLV


    def convert_to_cLorentzVector(self, p):
        """
          Convert four-vector of type LorentzVector to type cLorentzVector
          p four-vector [ak.Array]
          return: cLorentzVector(p)
        """
        print(" --- convert_to_cLorentzVector --- ")
        #vec = ak.zip(
        #    {
        #        "x": TComplex(p.x),
        #        "y": TComplex(p.y),
        #        "z": TComplex(p.z),
        #        "t": TComplex(p.t),
        #    },
        #    with_name="LorentzVector",
        #    behavior=vector.behavior,
        #)
        vec = np.zeros(4, dtype=TComplex)
        for mu in range(4):
            vec[mu] = TComplex(ak.to_numpy(self.get_component(p, mu)))
        #vec = ak.Array(vec)
        #vec = ak.from_regular(vec)
        #vec = vec[:,np.newaxis]
        print(f" >>> {vec}, {type(vec)}")
        return vec


    def convert_to_LorentzVector(self, p) -> ak.Array:
        """
          @brief Convert four-vector of type cLorentzVector to type LorentzVector. 
          Note: The imaginary part of the four-vector given as function argument is discarded in the conversion.
          @param p four-vector
          @return LorentzVector(p)
        """
        retVal = ak.zip(
            {
                "x": p[0].Re() if not isinstance(p[0], np.ndarray) else ak.Array(p[0]),
                "y": p[1].Re() if not isinstance(p[1], np.ndarray) else ak.Array(p[1]),
                "z": p[2].Re() if not isinstance(p[2], np.ndarray) else ak.Array(p[2]),
                "t": p[3].Re() if not isinstance(p[3], np.ndarray) else ak.Array(p[3]),
            },
            with_name = "LorentzVector",
            behavior = coffea.nanoevents.methods.vector.behavior,
        )
        return retVal


    
    def star(self, p) -> np.array:
        """
          Components must be awkward array.
          Numpy array is not allowed here.
          returns:
            awkward array of Lorentz vectors [each components are ak.array of complex conjugate objects]
        """
        #out = ak.zip(
        #    {
        #        "x": p.x.Conjugate(),
        #        "y": p.y.Conjugate(),
        #        "z": p.z.Conjugate(),
        #        "t": P.t.Conjugate()
        #    },
        #    with_name="LorentzVector",
        #    behavior=vector.behavior
        #)
        #return out
        retVal = np.zeros(4, dtype=TComplex)
        for mu in range(4):
            retVal[mu] = p[mu].Conjugate()
        return retVal

        

    def get_component(self, p, mu) -> ak.Array:
        """
          @brief Return certain component of given four-vector
          @param p given four-vector
          @param mu component to be returned [1]
          
          [1] use 0 for px
              use 1 for py
              use 2 for pz
              use 3 for energy
          
          @return p[mu]
        """
        if   mu == 0:    return p.x
        elif mu == 1:    return p.y
        elif mu == 2:    return p.z
        elif mu == 3:    return p.t
        else: assert False


    
    def comp_Pi(self, J, N):
        """
          J: numpy array of complex LorentzVector
             np.array([cLVx, cLVy, cLVz, cLVt])
             cLVx -> TComplex obejct -> each object Re / Im: np.array([[],[],[],....,[]])
          N: awkward array of LorentzVector
             ak.Array([{x:,y:,z:,t:},....]) 
        """
        print("  --- comp_Pi --- ")
        # LV -> complex LV
        # Now same as J
        cN = self.convert_to_cLorentzVector(N) 
        Jstar = self.star(J)

        JstarXN = Jstar.dot(self.g@cN)
        JXN     = J.dot(self.g@cN)
        JstarXJ = Jstar.dot(self.g@J)

        #print(f" >>> JstarXN: {JstarXN.Re()}, {JstarXN.Im()}")
        #print(f" >>> JXN: {JXN.Re()}, {JXN.Im()}")
        #print(f" >>> JstarXJ: {JstarXJ.Re()}, {JstarXJ.Im()}")

        cLV_1 = self.mult_arr_cLV(JstarXN,J)
        #print(f"cLV_1: {cLV_1}")
        cLV_2 = self.mult_arr_cLV(JXN,Jstar)
        #print(f"cLV_2: {cLV_2}")
        cLV_3 = self.mult_arr_cLV(JstarXJ,cN)
        #print(f"cLV_3: {cLV_3}")
        
        cLV = (self.mult_arr_cLV(JstarXN,J) + self.mult_arr_cLV(JXN,Jstar) - self.mult_arr_cLV(JstarXJ,cN))*2.0

        print(f" >>> cLV: {cLV}")
        
        #retVal = self.convert_to_LorentzVector(2.*(JstarXN*J + JXN*Jstar - JstarXJ*cN))
        retVal = self.convert_to_LorentzVector(cLV)
        print(f" >>> retVal: {retVal}")
        return retVal


    def mult_arr_cLV(self, arr, cLV):
        retVal = np.zeros(4, dtype=TComplex)
        for i in range(4):
            retVal[i] = arr * cLV[i]
        return retVal
            
    
    def comp_Pi5(self, J, N):
        """
          J: numpy array of complex LorentzVector
             np.array([cLVx, cLVy, cLVz, cLVt])
             cLVx -> TComplex obejct -> each object Re / Im: np.array([[],[],[],....,[]])
          N: awkward array of LorentzVector
             ak.Array([{x:,y:,z:,t:},....]) 
        """
        print("  --- comp_Pi 5 --- ")
        cN = self.convert_to_cLorentzVector(N)
        Jstar = self.star(J)

        #J_     = [J[0], J[1], J[2], J[3]]
        #Jstar_ = [Jstar[0], Jstar[1], Jstar[2], Jstar[3]]
        #cN_    = [cN[0], cN[1], cN[2], cN[3]]
        
        ##nev = ak.to_numpy(J.x).shape[0]
        #vProd = np.zeros(4, dtype=TComplex)
        ##vProd = self.get_epsilon(0,0,0,0)*Jstar[0]*J[0]*cN[0]
        #for mu in range(4):
        #    for nu in range(4):
        #        for rho in range(4):
        #            for sigma in range(4):
        #                epsilon = self.get_epsilon(mu, nu, rho, sigma)
        #                vProd[mu] += epsilon*Jstar[nu]*J[rho]*cN[sigma]

        vProd = np.zeros(4, dtype=TComplex)
        vProd[0] = self.get_epsilon(0,0,0,0)*Jstar[0]*J[0]*cN[0]
        vProd[1] = self.get_epsilon(1,0,0,0)*Jstar[0]*J[0]*cN[0]
        vProd[2] = self.get_epsilon(2,0,0,0)*Jstar[0]*J[0]*cN[0]
        vProd[3] = self.get_epsilon(3,0,0,0)*Jstar[0]*J[0]*cN[0]
        for nu in range(1,4):
            for rho in range(1,4):
                for sigma in range(1,4):
                    vProd[0] = vProd[0] + self.get_epsilon(0, nu, rho, sigma)*Jstar[nu]*J[rho]*cN[sigma]
                    vProd[1] = vProd[1] + self.get_epsilon(1, nu, rho, sigma)*Jstar[nu]*J[rho]*cN[sigma]
                    vProd[2] = vProd[2] + self.get_epsilon(2, nu, rho, sigma)*Jstar[nu]*J[rho]*cN[sigma]
                    vProd[3] = vProd[3] + self.get_epsilon(3, nu, rho, sigma)*Jstar[nu]*J[rho]*cN[sigma]
                        
        vProd = self.g@vProd
        #retVal = ak.zip(
        #    {
        #        "x": 2.*vProd[0].Im(),
        #        "y": 2.*vProd[1].Im(),
        #        "z": 2.*vProd[2].Im(),
        #        "t": 2.*vProd[3].Im()
        #    },
        #    with_name="LorentzVector",
        #    behavior=vector.behavior,            
        #)
        retVal = ak.zip(
            {
                "x": ak.Array(2.*vProd[0].Im()),
                "y": ak.Array(2.*vProd[1].Im()),
                "z": ak.Array(2.*vProd[2].Im()),
                "t": ak.Array(2.*vProd[3].Im())
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,            
        )
        print(f"retVal: {retVal}")
        return retVal


    def sgn(self, x):
        """
          @brief Return sign of integer value given as function argument
          @param x function argument
          @return sign(x)
        """
        if x > 0: return +1
        if x < 0: return -1
        return 0


    def get_epsilon(self, mu, nu, rho, sigma):
        """
          @brief Compute Levi-Civita symbol epsilon^{mu,nu,rho,sigma},
                 cf. https://en.wikipedia.org/wiki/Levi-Civita_symbol
                 (Section "Levi-Civita tensors")
          @params mu, nu,rho, sigma indices of Levi-Civita symbol
          @return epsilon^{mu,nu,rho,sigma}
          CV: formula for computation of four-dimensional Levi-Civita symbol taken from
               https://en.wikipedia.org/wiki/Levi-Civita_symbol
             (Section "Definition" -> "Generalization to n dimensions")
        """
        """
        a = 0*np.array[4]
        a[0] = mu
        a[1] = nu
        a[2] = rho
        a[3] = sigma
        epsilon = 1.
        for i in range(4):
            for j in range(4):
                epsilon *= self.sgn(a[j] - a[i])
        return epsilon
        """
        a = [mu, nu, rho, sigma]
        epsilon = 1.
        for i in range(4):
            for j in range(i+1,4):
                epsilon *= self.sgn(a[j] - a[i])
        return epsilon

    
    def kdash(self, si, mj, mk):
        """
         @brief Compute "decay momentum",
                given by bottom line of Eq. (A6) in Phys.Rev.D 61 (2000) 012002
         @param si (mass of two-pion system that forms resonance)^2
         @param mj mass of first  pion that forms resonance [1]
         @param mk mass of second pion that forms resonance [1]
        
         [1] the ordering of the two pions does not matter
        
         @return k'_i
        """
        retVal = np.sqrt((si - np.power(mj + mk, 2))*(si - np.power(mj - mk, 2)))/(2.0*np.sqrt(si))
        return retVal


    def Gamma(self, m0, Gamma0, si, mj, mk, L):
        """
        # @brief Compute "running width" of intermediate rho(770), rho(1450), f2(1270), sigma, and f0(1370)
        # resonances,
        #        given by bottom line of Eq. (A7) in Phys.Rev.D 61 (2000) 012002
        # @param m0 nominal mass of resonance
        # @param Gamma0 nominal width of resonance
        # @param si (mass of two-pion system that forms resonance)^2
        # @param mj mass of first  pion that forms resonance [1]
        # @param mk mass of second pion that forms resonance [1]
        # @param L angular momentum of resonance (s-wave = 0, p-wave=1, d-wave=2)
        #
        #  [1] the ordering of the two pions does not matter
        #
        # @return Gamma^{Y,L}(s_i)
        """
        kdashi = self.kdash(si, mj, mk)
        kdash0 = self.kdash(np.power(m0, 2), mj, mk)
        
        retVal = Gamma0*np.power(kdashi/kdash0, 2*L + 1)*m0/np.sqrt(si)
        return retVal


    def BreitWigner(self, m0, Gamma0, si, mj, mk, L) -> TComplex:
        """
          @brief Compute Breit-Wigner function of intermediate rho(770), rho(1450), f2(1270), sigma, and f0(1370) resonances,
                 given by top line of Eq. (A7) in Phys.Rev.D 61 (2000) 012002
          @param m0 nominal mass of resonance
          @param Gamma0 nominal width of resonance
          @param si (mass of two-pion system that forms resonance)^2
          @param mj mass of first  pion that forms resonance [1]
          @param mk mass of second pion that forms resonance [1]
          @param L angular momentum of resonance (s-wave = 0, p-wave=1, d-wave=2)
         
           [1] the ordering of the two pions does not matter
         
          @return B^{L}_{Y}(s_i)
        args:
          m0: scalar
          Gamma0: Scalar
          si: awkward array [mass]
          mj: scalar
          mk: scalar
          L: scalar
        return:
          Array of TComplex objects [awkward or numpy??? probably should be numpy :| ]
        """
        print(" --- BreitWigner --- ")
        num = -np.power(m0, 2)
        real = si + num
        imag = self.Gamma(m0, Gamma0, si, mj, mk, L)*m0
        #print(f" >>> Num  : {num}")
        #print(f" >>> Real : {real}")
        #print(f" >>> Imag : {imag}")
        ##denom = TComplex(real, imag)
        denom = TComplex(ak.to_numpy(real), ak.to_numpy(imag))
        ##denom = TComplex(ak.to_numpy(si - np.power(m0, 2)), ak.to_numpy(m0*self.Gamma(m0, Gamma0, si, mj, mk, L)))
        #print(f" >>> denom: {type(denom)}\t{denom.Re()},{denom.Re().shape}\t{denom.Im()},{denom.Im().shape}")
        retval = TComplex(num, 0.)/denom
        #print(f" >>> retval: {type(retval)}\t{retval.Re()},{retval.Re().shape}\t{retval.Im()},{retval.Im().shape}")

        return retval


    def m_a1(self, s: ak.Array) -> float:
        return self.m0_a1


    def Gamma_a1(self, s: ak.Array) -> ak.Array:
        print(f"s: {s}")
        
        mask1       = (s - self.Gamma_a1_vs_s[0][0]) <= 0
        #print(f"mask1: {mask1}")
        mask1result = self.Gamma_a1_vs_s[0][1]*ak.ones_like(s)
        #print(f"mask1result: {mask1result}")
        mask2       = (s - self.Gamma_a1_vs_s[-1][0]) >= 0
        #print(f"mask2: {mask2}")
        mask2result = self.Gamma_a1_vs_s[-1][1]*ak.ones_like(s) 
        #print(f"mask2result: {mask2result}")
        
        gamma_left  = ak.Array([a[0] for a in self.Gamma_a1_vs_s])
        #print(f"gamma_left: {gamma_left}")
        gamma_right = ak.Array([a[1] for a in self.Gamma_a1_vs_s])
        #print(f"gamma_right: {gamma_right}")
        
        gammaS, gammaL = ak.broadcast_arrays(ak.to_numpy(s), ak.to_numpy(gamma_left))
        _     , gammaR = ak.broadcast_arrays(ak.to_numpy(s), ak.to_numpy(gamma_right))

        gammaS = ak.from_regular(gammaS)
        gammaL = ak.from_regular(gammaL)
        gammaR = ak.from_regular(gammaR)
        print(gammaS, gammaL, gammaR)
        
        diffLS = gammaS - gammaL
        print(f"diffLS: {diffLS}")
        mask_low  = diffLS >= 0.0
        mask_high = diffLS <= 0.0
        print(f"mask_low: {mask_low}")
        print(f"mask_high: {mask_high}")
        
        idx_low    = ak.local_index(gammaL)[mask_low][:,-1:]
        idx_high   = ak.local_index(gammaL)[mask_high][:,0:1]
        print(f"idx_low: {idx_low}, {idx_low.type}")
        print(f"idx_high: {idx_high}, {idx_high.type}")

        
        #print(idx_low[:,-1:])
        #print(idx_high[:,0:1])
        print(gammaL[idx_low])
        print(gammaL[idx_high])
        
        #s_lo = gammaL[idx_low][:,-1][:,None]
        #s_hi = gammaL[idx_high][:,1][:,None]
        #Gamma_lo = gammaR[idx_low][:,-1][:,None]
        #Gamma_hi = gammaR[idx_high][:,1][:,None]
        s_lo = gammaL[idx_low]
        s_hi = ak.fill_none(ak.pad_none(gammaL[idx_high], 1, axis=-1), 0)
        Gamma_lo = gammaR[idx_low]
        Gamma_hi = ak.fill_none(ak.pad_none(gammaR[idx_high], 1, axis=-1), 0)

        #for i in range(35094):
        #    print(s_lo[i], s_hi[i], Gamma_lo[i], Gamma_hi[i])

        
        
        #print(ak.min(ak.num(s_lo, axis=1)),     ak.sum(ak.num(s_lo, axis=1)), s_lo[1816], s_lo[7363])
        #print(ak.min(ak.num(s_hi, axis=1)),     ak.sum(ak.num(s_hi, axis=1)), s_hi[1816], s_hi[7363])
        #print(ak.min(ak.num(Gamma_lo, axis=1)), ak.sum(ak.num(Gamma_lo, axis=1)), Gamma_lo[1816], Gamma_lo[7363])
        #print(ak.min(ak.num(Gamma_hi, axis=1)), ak.sum(ak.num(Gamma_hi, axis=1)), Gamma_hi[1816], Gamma_hi[7363])

        print(ak.argsort(ak.num(s_hi, axis=1), ascending=True))
        
        #print((ak.to_numpy(s)).shape, (ak.to_numpy(s_lo)).shape, (ak.to_numpy(s_hi)).shape)
        print(s.type, s_lo.type, s_hi.type)
        
        retVal = ak.where(mask1,
                          mask1result,
                          ak.where(mask2,
                                   mask2result,
                                   (Gamma_lo*(s_hi - s) + Gamma_hi*(s - s_lo)) / (s_hi - s_lo)
                                   )
                          )

        print(retVal)
        return retVal


    def BreitWigner_a1(self, s: ak.Array) -> TComplex:
        """
          return: awkward array of TComplex objects
        """
        m = self.m_a1(s)
        Gamma = self.Gamma_a1(s)

        num = -np.power(m, 2)
        denom = TComplex(s - np.power(m, 2), self.m0_a1*Gamma)

        retVal = TComplex(num, 0.)/denom
        return retVal
