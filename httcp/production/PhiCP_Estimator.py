import os
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

from httcp.production.PolarimetricA1 import PolarimetricA1


def GetPhiCP(
        p4hcandinfodict: dict, 
        method_leg1: str, 
        method_leg2: str,
        mode_leg1: str, 
        mode_leg2: str,
        **kwars
) -> ak.Array:
    #print(f"{method_leg1}-{method_leg2}")
    #print(f"{mode_leg1}-{mode_leg2}")
    """
    Inputs:
      p4hcandinfodict : full hcand dict
        e.g. 
          hcand : [h1 (i.e. e/mu/tau), h2 (i.e. tau)]
          p4hcandinfodict = {
            "p4h1"    : p4 of e/mu/tau i.e. the 1st component of hcand,
            "p4h1pi"  : p4 of charged pion if tau else empty,
            "p4h1pi0" : p4 of reconstructed pi0 if tau else empty,
            "p4h2"    : p4 of tau (always) i.e. the 2nd component of hcand,
            "p4h2pi"  : p4 of charged pion from tau,
            "p4h2pi0" : p4 of reconstructed pi0 from tau,
          }
      method_leg1     : "DP/PV/IP"
      method_leg2     : "DP/PV/IP"
      mode_leg1       : "e/mu/pi/rho/a1"
      mode_leg2       : "pi/rho/a1" because the 2nd component of hcand is always tauh
    Output:
      Final PhiCP array for any methods / configurations
    Steps:
      Prepares the input vectors for PhiCP
      Compute the PhiCP
    """
    #from IPython import embed; embed()
    phicp_input_vec_dict = PrepareVecsForPhiCP(p4hcandinfodict,
                                               method_leg1,
                                               method_leg2, 
                                               mode_leg1,
                                               mode_leg2)
    phicp = ak.values_astype(ComputeAcopAngle(phicp_input_vec_dict), np.float32)
    phicp = ak.enforce_type(phicp, "var * float32")
    
    return phicp


def PrepareVecsForPhiCP(
        p4hcandinfodict: dict, 
        method_leg1: str, 
        method_leg2: str, 
        mode_leg1: str, 
        mode_leg2: str
) -> dict:
    """
    Inputs:
      p4hcandinfodict : full hcand dict
      method_leg1     : "DP/PV/IP"
      method_leg2     : "DP/PV/IP"
      mode_leg1       : "e/mu/pi/rho/a1"
      mode_leg2       : "pi/rho/a1"
    Output:
      P1, R1, P2, R2, pi_phase_shift =  y1*y2
      i = 1,2
      Impact Parameter (IP) method:
        Pi = momemtum of tau decau product (e/mu/pi) at zero momemtum frame (ZMF)
        Ri = IP vector at ZMF
        yi = phase shift = 0 for this method 

      Decay Plane (DP) or neutral pion method:
        Pi = momentum of charged pion at ZMF
        Ri = momentum of the neutral pion at ZMF
        yi = (Epi - Epi0)/(Epi + Epi0) 

      Polarimetric Vector (PV) method:
        Pi = tau momentum at ZMF
        Ri = polarimetric vector at ZMF
        yi = 0 for this method
    Steps:
      Prepare input vectors for PhiCP calculation
    """
    p4h1    = p4hcandinfodict["p4h1"]
    p4h1pi  = p4hcandinfodict["p4h1pi"]
    p4h1pi0 = p4hcandinfodict["p4h1pi0"]
    p4h2    = p4hcandinfodict["p4h2"]
    p4h2pi  = p4hcandinfodict["p4h2pi"]
    p4h2pi0 = p4hcandinfodict["p4h2pi0"]
    
    _P1, _R1, _y1, _c1 = _prepareVecs(p4h1, p4h1pi, p4h1pi0, method_leg1, mode_leg1)
    _P2, _R2, _y2, _c2 = _prepareVecs(p4h2, p4h2pi, p4h2pi0, method_leg2, mode_leg2)

    # Get the boost vector
    boostv = _getBoost(_P1, _P2)

    # Boost the vectors
    if method_leg1 == "PV":
      P1, R1, y1 = _boostVec_and_PV(boostv, _P1, _y1, method_leg1, mode_leg1, p4h1, p4h1pi, p4h1pi0)
    else:
      P1, R1, y1 = _boostVecs(boostv, _P1, _R1, _y1, method_leg1, mode_leg1)

    if method_leg2 == "PV":
      P2, R2, y2 = _boostVec_and_PV(boostv, _P2, _y2, method_leg2, mode_leg2, p4h2, p4h2pi, p4h2pi0)
    else:
      P2, R2, y2 = _boostVecs(boostv, _P2, _R2, _y2, method_leg2, mode_leg2)

    return {"P1" : P1, "R1" : R1, 
            "P2" : P2, "R2" : R2,
            "Y1" : y1, "Y2" : y2,
            "C1" : _c1, "C2": _c2}


def _boostVec_and_PV(
        boostv: ak.Array,
        _P: ak.Array,
        _y: ak.Array,
        method_leg: str,
        mode_leg: str,
        p4_hcand: ak.Array,
        p4_hcand_pi: ak.Array,
        p4_hcand_pi0: ak.Array
) -> tuple[ak.Array, ak.Array, ak.Array]:
    """
    Boost the vectors and calculate the polarimetric vector
    """
    P = _P.boost(boostv.negative())
    R = None
    Y = _y
    if mode_leg == "pi":
        R = p4_hcand_pi.boost(boostv.negative()).pvec
    elif mode_leg == "rho":
        pi  = p4_hcand_pi.boost(boostv.negative())
        pi0 = p4_hcand_pi0.boost(boostv.negative())
        q   = pi.subtract(pi0)
        N   = P.subtract(pi.add(pi0))
        pv  = (((2*(q.dot(N))*q.pvec).subtract(q.mass2*N.pvec)))
        R = pv #.unit
    elif mode_leg == "a1":
        os_pi_HRF   = p4_hcand_pi[:, 0:1].boost(boostv.negative())
        ss1_pi_HRF  = p4_hcand_pi[:, 1:2].boost(boostv.negative()) 
        ss2_pi_HRF  = p4_hcand_pi[:, 2:3].boost(boostv.negative()) 
        a1pol       = PolarimetricA1(P,
                                     os_pi_HRF,
                                     ss1_pi_HRF,
                                     ss2_pi_HRF,
                                     p4_hcand.charge)
        R = -a1pol.PVC().pvec #.pvec #.unit
    else:
        raise RuntimeError(f"Wrong mode: {mode_leg}")

    return P, R, Y

def _boostVecs(
        boostv: ak.Array,
        _P: ak.Array,
        _R: ak.Array,
        _y: ak.Array,
        method_leg: str, 
        mode_leg: str, 
) -> tuple[ak.Array, ak.Array, ak.Array]:
    """
    Boost the vectors
    """
    P = _P.boost(boostv.negative())
    R = _R.boost(boostv.negative())
    Y = _y

    return P.pvec, R.pvec, Y
    #return P, R, Y


def ComputeAcopAngle(vecsdict):
    """
    Geometrical estimation of PhiCP
    For PVPV:
       P1/P2  --> full tau p4
       R1/R2  --> k i.e. full tau x polarimetric vector
       H1/H2  --> h i.e. polarimetric vector
    For DPDP: 
       P1/P2  --> vecPiPlus
       R1/R2  --> vecPiZeroPlustransv
       H1/H2  --> R1/R2
       Y1     --> phase-shift for leg 1
       Y2     --> phase-shift for leg 2
    """
    #embed()
    # assert len(vecsdict) == 10, "input dict to ComputeAcopAngle does not have proper structure"
    P1 = vecsdict["P1"]
    R1 = vecsdict["R1"]
    P2 = vecsdict["P2"]
    R2 = vecsdict["R2"]
    Y1 = vecsdict["Y1"]
    Y2 = vecsdict["Y2"]
    C1 = vecsdict["C1"]
    C2 = vecsdict["C2"]
    
    # Calculate the phase shift
    Y = Y1*Y2

    # Calculate R_perp
    R1_perp = R1 - R1.dot(P1.unit)*P1.unit
    R2_perp = R2 - R2.dot(P2.unit)*P2.unit
    """
    O_sign = (R2_perp.cross(R1_perp)).dot(P1)
    """
    acop = np.arccos(R1_perp.unit.dot(R2_perp.unit))
    
    # to get the sign
    Pm = ak.where(C1 < 0, P1, P2)
    Rm_perp = ak.where(C1 < 0, R1_perp, R2_perp)
    Pp = ak.where(C1 < 0, P2, P1)
    Rp_perp = ak.where(C1 < 0, R2_perp, R1_perp)
    O_sign = (Rp_perp.cross(Rm_perp)).dot(Pm)
    
    # Apply the shift if sign is negative
    acop  = ak.where(O_sign < 0.0, 2.*np.pi - acop, acop)

    acop  = ak.where(Y < 0.0, acop + np.pi, acop) # effective for DP method only
    
    #Map  angles into [0,2pi] interval
    acop = ak.where(acop > 2.*np.pi, acop - 2.* np.pi, acop) 
    acop = ak.where(acop < 0,        acop + 2.* np.pi, acop)

    return acop


def _getBoost(
      P1: ak.Array,
      P2: ak.Array,
)-> ak.Array:
  """
  A function to get the boostvec of the corresponding frame
  """
  frame = P1 + P2

  return frame.boostvec


def _prepareVecs(
        p4_hcand: ak.Array, # tau momentum
        p4_hcand_pi: ak.Array, # pion momentum
        p4_hcand_pi0: ak.Array, 
        leg_method: str, 
        leg_mode: str
) -> tuple[ak.Array, ak.Array, ak.Array]:
    """
    A private function to be used in PrepareVecsForPhiCP
    This mainly returns the required vectors for different modes but not boosted
    Returned vectors will further be modified in the reStructure step
    Output:
    _P, R, _y
    Impact Parameter (IP) method:
      Pi = momemtum of tau decau product (e/mu/pi) 
      Ri = IP vector 
      yi = phase shift = 0 for this method 

    Decay Plane (DP) or neutral pion method:
      Pi = momentum of charged pion 
      Ri = momentum of the neutral pion 
      yi = (Epi - Epi0)/(Epi + Epi0) 

    Polarimetric Vector (PV) method:
      Pi = tau momentum 
      Ri = polarimetric vector 
      yi = 0 for this method
      
    """
    _P = None
    _R = None
    _y = None

    if leg_method == "IP":
      if leg_mode == "e" or leg_mode == "mu":
          _P = p4_hcand # e/mu the tau is already associated to the e/mu
          _R = _calculateIP(p4_hcand)
          _y = ak.ones_like(p4_hcand.energy)
      elif leg_mode == "pi":
          _P = p4_hcand_pi
          _R = _calculateIP(p4_hcand)
          _y = ak.ones_like(p4_hcand.energy)
      else:
          raise RuntimeError(f"Wrong mode : {leg_mode}")
   
    elif leg_method == "DP":
      # https://github.com/alebihan/IPHCProductionTools/blob/master/HiggsCPinTauDecays/TauDecaysInterface/src/SCalculator.cc#L363-L366
      if leg_mode == "rho":
          _P = p4_hcand_pi
          _R = p4_hcand_pi0
          #_y = (p4_hcand_pi.energy - p4_hcand_pi0.energy)/(p4_hcand_pi.energy + p4_hcand_pi0.energy)
          _y = _getshift(_P, _R)
      elif leg_mode == "a1":
          _P = _get_pi_a1_DP(p4_hcand_pi)
          _R = p4_hcand_pi[:,:1]
          #_y = (p4_hcand_pi[:,:1].energy - p4_hcand_pi[:,1:2].energy)/(p4_hcand_pi[:,:1].energy + p4_hcand_pi[:,1:2].energy)
          _y = _getshift(_P, _R) #(_P.energy - _R.energy)/(_P.energy + _R.energy)
      else:
          raise RuntimeError(f"Wrong mode : {leg_mode}")
    
    elif leg_method == "PV":
      _P = p4_hcand
      _R = None # polarimetic vector need to be calculated with boost
      _y = ak.ones_like(p4_hcand.energy)  # ak aray of 1.0
    else:
        raise RuntimeError(f"Wrong {leg_method}")

    return _P, _R, _y, p4_hcand.charge


def _calculateIP(hcand: ak.Array) -> ak.Array:
    """
    Calculate the impact parameter vector and retrun it
    """
    IP = ak.zip(
        {
            "x":hcand.IPx,
            "y":hcand.IPy,
            "z":hcand.IPz,
            "t":ak.zeros_like(hcand.IPz)
        },
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior)
    return IP


def _getshift(p4P, p4R):
    return (p4P.energy - p4R.energy)/(p4P.energy + p4R.energy)


def _get_pi_a1_DP(p4_pi):
    Minv1 = (p4_pi[:,:1] + p4_pi[:,1:2]).mass
    Minv2 = (p4_pi[:,:1] + p4_pi[:,2:3]).mass
    Pi = ak.where(np.abs(0.77526-Minv1) < np.abs(0.77526-Minv2), p4_pi[:,1:2], p4_pi[:,2:3])
    return Pi
