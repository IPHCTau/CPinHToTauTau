from columnflow.util import maybe_import

ak = maybe_import("awkward")

from rich.console import Console
from rich.tree import Tree


def leptau_category_flow(tag, events):
    flowdict = {}
    
    flowdict["tau is real"] = events.is_real_2
    flowdict["tau is fake"] = events.is_fake_2
    
    flowdict["A  : is_ss + is_iso_2    + is_b_veto + is_low_mt "]  = ~events.is_os & events.is_iso_2 & events.is_b_veto & events.is_low_mt
    flowdict["B  : is_ss + is_noniso_2 + is_b_veto + is_low_mt "]  = ~events.is_os & ~events.is_iso_2 & events.is_b_veto & events.is_low_mt
    flowdict["A0 : is_os + is_iso_2    + has_b     + is_low_mt "]  = events.is_os & events.is_iso_2 & ~events.is_b_veto & events.is_low_mt
    flowdict["B0 : is_os + is_noniso_2 + has_b     + is_low_mt "]  = events.is_os & ~events.is_iso_2 & ~events.is_b_veto & events.is_low_mt
    flowdict["A1 : is_os + is_iso_2    + is_b_veto + is_high_mt"]  = events.is_os & events.is_iso_2 & events.is_b_veto & ~events.is_low_mt
    flowdict["B1 : is_os + is_noniso_2 + is_b_veto + is_high_mt"]  = events.is_os & ~events.is_iso_2 & events.is_b_veto & ~events.is_low_mt
    flowdict["D  : is_os + is_iso_2    + is_b_veto + is_low_mt "]  = events.is_os & events.is_iso_2 & events.is_b_veto & events.is_low_mt
    flowdict["C  : is_os + is_noniso_2 + is_b_veto + is_low_mt "]  = events.is_os & ~events.is_iso_2 & events.is_b_veto & events.is_low_mt
    
    flowdict["DM = 0 "]  = events.is_pi_2
    flowdict["DM = 1 "]  = events.is_rho_2
    flowdict["DM = 2 "]  = events.is_a1_1pr_2pi0_2
    flowdict["DM = 10"]  = events.is_a1_3pr_0pi0_2
    flowdict["DM = 11"]  = events.is_a1_3pr_1pi0_2
    
    
    root = Tree(f"[yellow]{tag} : {ak.count(events.event)}")
    
    for a in ["tau is real", "tau is fake"]:
        _temp = root.add(f"[cyan]{a} : {ak.sum(flowdict[a])}")
        for b in ["A  : is_ss + is_iso_2    + is_b_veto + is_low_mt ",
                  "B  : is_ss + is_noniso_2 + is_b_veto + is_low_mt ",
                  "A0 : is_os + is_iso_2    + has_b     + is_low_mt ",
                  "B0 : is_os + is_noniso_2 + has_b     + is_low_mt ",
                  "A1 : is_os + is_iso_2    + is_b_veto + is_high_mt",
                  "B1 : is_os + is_noniso_2 + is_b_veto + is_high_mt",
                  "D  : is_os + is_iso_2    + is_b_veto + is_low_mt ",
                  "C  : is_os + is_noniso_2 + is_b_veto + is_low_mt "]:
            _temp2 = _temp.add(f"[green]{b} : {ak.sum(flowdict[a] & flowdict[b])}")
            for c in ["DM = 0 ", "DM = 1 ", "DM = 2 ", "DM = 10", "DM = 11"]:
                _temp2.add(f"[bright_white]{c} : {ak.sum(flowdict[a] & flowdict[b] & flowdict[c])}")
    return root


def tautau_category_flow(tag, events):
    flowdict = {}
    
    flowdict["leading tau is real"] = events.is_real_1
    flowdict["leading tau is fake"] = events.is_fake_1
    
    flowdict["A  : is_ss + is_iso_1    + is_iso_2    "] = ~events.is_os & events.is_iso_1 & events.is_iso_2
    flowdict["B  : is_ss + is_noniso_1 + is_iso_2    "] = ~events.is_os & ~events.is_iso_1 & events.is_iso_2
    flowdict["A0 : is_ss + is_noniso_1 + is_noniso_2 "] = ~events.is_os & events.is_iso_1 & ~events.is_iso_2
    flowdict["B0 : is_ss + is_noniso_1 + is_noniso_2 "] = ~events.is_os & ~events.is_iso_1 & ~events.is_iso_2
    flowdict["D0 : is_os + is_iso_1    + is_noniso_2 "] = events.is_os & events.is_iso_1 & ~events.is_iso_2
    flowdict["C0 : is_os + is_noniso_1 + is_noniso_2 "] = events.is_os & ~events.is_iso_1 & ~events.is_iso_2
    flowdict["D  : is_os + is_iso_1    + is_iso_2    "] = events.is_os & events.is_iso_1 & events.is_iso_2
    flowdict["C  : is_os + is_noniso_1 + is_iso_2    "] = events.is_os & ~events.is_iso_1 & events.is_iso_2
    
    flowdict["DM1 = 0"]     = events.is_pi_1
    flowdict["DM1 = 1"]     = events.is_rho_1
    flowdict["DM1 = 2"]     = events.is_a1_1pr_2pi0_1
    flowdict["DM1 = 10"]    = events.is_a1_3pr_0pi0_1
    flowdict["DM1 = 11"]    = events.is_a1_3pr_1pi0_1
    #flowdict["DMs = 0-0"]   = events.is_pi_1 & events.is_pi_2
    #flowdict["DMs = 0-1"]   = ((events.is_pi_1 & events.is_rho_2) | (events.is_pi_2 & events.is_rho_1))
    #flowdict["DMs = 0-2"]   = ((events.is_pi_1 & events.is_a1_1pr_2pi0_2) | (events.is_pi_2 & events.is_a1_1pr_2pi0_1))
    #flowdict["DMs = 0-10"]  = ((events.is_pi_1 & events.is_a1_3pr_0pi0_2) | (events.is_pi_2 & events.is_a1_3pr_0pi0_1))
    #flowdict["DMs = 0-11"]  = ((events.is_pi_1 & events.is_a1_3pr_1pi0_2) | (events.is_pi_2 & events.is_a1_3pr_1pi0_1))
    #flowdict["DMs = 1-1"]   = events.is_rho_1 & events.is_rho_2
    #flowdict["DMs = 1-2"]   = ((events.is_rho_1 & events.is_a1_1pr_2pi0_2) | (events.is_rho_2 & events.is_a1_1pr_2pi0_1))
    #flowdict["DMs = 1-10"]  = ((events.is_rho_1 & events.is_a1_3pr_0pi0_2) | (events.is_rho_2 & events.is_a1_3pr_0pi0_1))
    #flowdict["DMs = 1-11"]  = ((events.is_rho_1 & events.is_a1_3pr_1pi0_2) | (events.is_rho_2 & events.is_a1_3pr_1pi0_1))
    #flowdict["DMs = 2-2"]   = events.is_a1_1pr_2pi0_1 & events.is_a1_1pr_2pi0_2
    #flowdict["DMs = 2-10"]  = ((events.is_a1_1pr_2pi0_1 & events.is_a1_3pr_0pi0_2) | (events.is_a1_1pr_2pi0_2 & events.is_a1_3pr_0pi0_1))
    #flowdict["DMs = 2-11"]  = ((events.is_a1_1pr_2pi0_1 & events.is_a1_3pr_1pi0_2) | (events.is_a1_1pr_2pi0_2 & events.is_a1_3pr_1pi0_1))
    #flowdict["DMs = 10-10"] = events.is_a1_3pr_0pi0_1 & events.is_a1_3pr_0pi0_2
    #flowdict["DMs = 10-11"] = ((events.is_a1_3pr_0pi0_1 & events.is_a1_3pr_1pi0_2) | (events.is_a1_3pr_0pi0_2 & events.is_a1_3pr_1pi0_1))
    #flowdict["DMs = 11-11"] = events.is_a1_3pr_1pi0_1 & events.is_a1_3pr_1pi0_2
    
    root = Tree(f"[yellow]{tag} : {ak.count(events.event)}")
    
    for a in ["leading tau is real", "leading tau is fake"]:
        _temp = root.add(f"[cyan]{a} : {ak.sum(flowdict[a])}")
        for b in ["A  : is_ss + is_iso_1    + is_iso_2    ",
                  "B  : is_ss + is_noniso_1 + is_iso_2    ",
                  "A0 : is_ss + is_noniso_1 + is_noniso_2 ",
                  "B0 : is_ss + is_noniso_1 + is_noniso_2 ",
                  "C0 : is_os + is_noniso_1 + is_noniso_2 ",
                  "D0 : is_os + is_iso_1    + is_noniso_2 ",
                  "C  : is_os + is_noniso_1 + is_iso_2    ",
                  "D  : is_os + is_iso_1    + is_iso_2    ",
                  ]:
            _temp2 = _temp.add(f"[green]{b} : {ak.sum(flowdict[a] & flowdict[b])}")
            for c in [
                    "DM1 = 0", "DM1 = 1", "DM1 = 2", "DM1 = 10", "DM1 = 11",
                    #"DMs = 0-0", "DMs = 0-1", "DMs = 0-2", "DMs = 0-10", "DMs = 0-11",
                    #"DMs = 1-1", "DMs = 1-2", "DMs = 1-10", "DMs = 1-11",
                    #"DMs = 2-2", "DMs = 2-10", "DMs = 2-11",
                    #"DMs = 10-10", "DMs = 10-11",
                    #"DMs = 11-11"
            ]:
                _temp2.add(f"[bright_white]{c} : {ak.sum(flowdict[a] & flowdict[b] & flowdict[c])}")
    return root




def category_flow(tag, events):
    # Create a console object
    console = Console(record=True)
    if tag == "etau" or tag == "mutau":
        root = leptau_category_flow(tag, events)
    elif tag == "tautau":
        root = tautau_category_flow(tag, events)
        
    console.print(root)
