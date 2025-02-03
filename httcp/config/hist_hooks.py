# coding: utf-8

"""
Histogram hooks.
"""

from __future__ import annotations

from collections import defaultdict

import os
import law
import order as od
import scinum as sn
import json
import pickle
import numpy as np
from collections import defaultdict

from columnflow.util import maybe_import, DotDict

np = maybe_import("numpy")
hist = maybe_import("hist")



logger = law.logger.get_logger(__name__)


def add_hist_hooks(config: od.Config) -> None:
    """
    Add histogram hooks to a configuration.
    """
    # helper to convert a histogram to a number object containing bin values and uncertainties
    # from variances stored in an array of values
    def hist_to_num(h: hist.Histogram, unc_name=str(sn.DEFAULT)) -> sn.Number:
        return sn.Number(h.values(), {unc_name: h.variances()**0.5})

    # helper to integrate values stored in an array based number object
    def integrate_num(num: sn.Number, axis=None) -> sn.Number:
        return sn.Number(
            nominal=num.nominal.sum(axis=axis),
            uncertainties={
                unc_name: (
                    (unc_values_up**2).sum(axis=axis)**0.5,
                    (unc_values_down**2).sum(axis=axis)**0.5,
                )
                for unc_name, (unc_values_up, unc_values_down) in num.uncertainties.items()
            },
        )

    ######################################################
    ###         PRODUCE FAKE FACTOR HISTOGRAMS         ###
    ######################################################
    def produce_fake_factor(task, hists):
        '''
        Produce fake factor histograms for ABCD or A0B0C0D0 categories
        FF = (data-mc)A/(data-mc)B
        FF0 = (data-mc)A0/(data-mc)B0
        The resulting fake factors histograms are stroe in a pickle file.

        *** N.B. 
        * This will give : TypeError: 'NoneType' object is not iterable : Error at the end
        as this finction does not return hist. But the FF pkl files will be saved. 
        * One can use one histogram only, e.g. hcand_1_pt
        '''
        # Define if we are calculating the fake factor for A0B0C0D0 or ABCD categories
        iszero_category = False # True for A0B0C0D0 categories

        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd process, it will be used to store the fake factor histograms
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            print("no qcd process")
            return hists


        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups: A, B, C, D of A0, B0, C0, D0 for each DM and Njet category
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)


        dms = ["tau1a1DM11", "tau1a1DM10", "tau1a1DM2", "tau1pi", "tau1rho"]  # Decay modes
        njets = ["has0j", "has1j", "has2j"]  # Jet multiplicity

        # Loop over all categories and create a QCD group for each DM and Njet category
        for dm in dms:
            for njet in njets:
                for cat_id in category_ids:
                    cat_inst = config.get_category(cat_id)

                    if iszero_category: # A0B0C0D0 categories
                        if cat_inst.has_tag({"ss", "iso1", "noniso2", njet, dm}, mode=all): # cat A0
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_iso = cat_inst
                        elif cat_inst.has_tag({"ss", "noniso1", "noniso2", njet, dm}, mode=all): # cat BB
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_noniso = cat_inst
                    else: # ABCD categories
                        if cat_inst.has_tag({"ss", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat A 
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_iso = cat_inst
                        elif cat_inst.has_tag({"ss", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat B
                            qcd_groups[f"dm_{dm}_njet_{njet}"].ss_noniso = cat_inst

        # Get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 2]

        # Nothing to do if there are no complete groups, you need A and B to estimate FF
        if not complete_groups:
            print("no complete groups")
            return hists
        
        # Sum up mc and data histograms, stop early when empty, this is done for all categories
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy()) # sum all MC histograms, the hist object here contains all categories
        data_hist = sum(data_hists[1:], data_hists[0].copy()) # sum all data histograms, the hist object here contains all categories

        for gidx, group_name in enumerate(complete_groups):
            group = qcd_groups[group_name]
            # Get the corresponding histograms of the id, if not present, create a zeroed histogram
            get_hist = lambda h, region_name: (
                h[{"category": hist.loc(group[region_name].id)}]
                if group[region_name].id in h.axes["category"]
                else hist.Hist(*[axis for axis in (h[{"category": [0]}] * 0).axes if axis.name != 'category'])
            )

            # Get the corresponding histograms and convert them to number objects,
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc") # MC in region A
            ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso"), "ss_iso_data") # Data in region A
            ss_noniso_mc  = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc") # MC in region B
            ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso"), "ss_noniso_data") # Data in region B


            # calculate the pt-dependent fake factor
            fake_factor = ((ss_iso_data - ss_iso_mc) / (ss_noniso_data - ss_noniso_mc))[:, None] # FF = (data-mc)A/(data-mc)B
            fake_factor_values = np.squeeze(np.nan_to_num(fake_factor()), axis=0) # get the values of the fake factor
            fake_factor_variances = fake_factor(sn.UP, sn.ALL, unc=True)**2 # get the uncertainties of the fake factor

            fake_factor_variances = fake_factor_variances[0]

            # Guaranty positive values of fake_factor
            neg_int_mask = fake_factor_values <= 0
            fake_factor_values[neg_int_mask] = 1e-5
            fake_factor_variances[neg_int_mask] = 0
        

            # create a hist clone of the data_hist
            ratio_hist = data_hist.copy()

            # fill the ratio histogram with the fake factor values in the first category
            ratio_hist.view().value[0, ...] = fake_factor_values
            ratio_hist.view().variance[0, ...] = fake_factor_variances

            # Save the fake factor histogram in a pickle file
            #path = "/eos/user/o/oponcet2/analysis/CP_dev/analysis_httcp/cf.PlotVariables1D/FF"
            path = f"{law.wlcg.WLCGFileSystem().base[0].split('root://eosuser.cern.ch')[-1]}/analysis_httcp/cf.PlotVariables1D/{config.campaign.name}/FF"
            hname = ratio_hist.axes[2].name
            
            # Ensure the folder exists
            if not os.path.exists(path):
                os.makedirs(path)
            with open(f"{path}/fake_factors_{hname}_{group_name}.pkl", "wb") as f:
                pickle.dump(ratio_hist, f)

 
    ######################################################
    ###           EXTRAPOLATE FAKE PROCESS             ###
    ######################################################
    def extrapolate_fake(task, hists):
        '''
        This is a multi task function that can create extrapolate the fake process of one region to another.
        - FF x B -> A : type_extrapolation = "AB"; also create ratio plot of DATA/MC of A region for closure correction
        - FF x C -> D : type_extrapolation = "CD" it's control plots
        - FF0 x C0 -> D0 : type_extrapolation = "C0D0"; also create ratio plot of DATA/MC oof D0 region for extrapolation correction
        '''
        # Choose the type of extrapolation
        #type_extrapolation = "CD" # "AB" or "CD" or "C0D0"
        type_extrapolation = config.x.regions_to_extrapolate_fake

        
        # Check if histograms are available
        if not hists:
            print("no hists")
            return hists

        # Get the qcd proces, this will be used as the fake process
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            print("no fake") 
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        ### Create a QCD group ofr each DM and Njet category
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)

        dms = ["tau1a1DM11", "tau1a1DM10", "tau1a1DM2", "tau1pi", "tau1rho"]  # Decay modes
        njets = ["has0j", "has1j", "has2j"]  # Jet multiplicity


        # Loop over all categories and create a QCD group for each DM and Njet category
        for dm in dms:
            for njet in njets:
                for cat_id in category_ids:
                    cat_inst = config.get_category(cat_id)

                    # CASE OF AB CLOSURE PLOTS
                    if type_extrapolation == "AB":
                        if cat_inst.has_tag({"ss", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat A 
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                        elif cat_inst.has_tag({"ss", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat B
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst

                    # CASE OF CD CONTROL PLOTS
                    elif type_extrapolation == "CD":
                        if cat_inst.has_tag({"os", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat D
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                        elif cat_inst.has_tag({"os", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat C
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst
                    
                    # CASE OF C0D0 CONTROL PLOTS
                    elif type_extrapolation == "C0D0":
                        if cat_inst.has_tag({"os", "iso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat D0
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_iso = cat_inst
                        elif cat_inst.has_tag({"os", "noniso1", njet, dm}, mode=all) and not cat_inst.has_tag("noniso2"): # cat C0
                            qcd_groups[f"dm_{dm}_njet_{njet}"].os_noniso = cat_inst   

        # Get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 2]
    
        # Nothing to do if there are no complete groups, you need C to apply Fake to D 
        if not complete_groups:
            print("no complete groups")
            return hists

        # Sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())
        
        # Start by copying the data hist and reset it, then fill it at specific category slices
        hists[qcd_proc] = qcd_hist = data_hist.copy().reset()
        mc_hist_incl = mc_hist.copy().reset()
        data_hist_incl = data_hist.copy().reset()
        os_iso_mc_incl = None
        os_iso_data_incl = None
        path = f"{law.wlcg.WLCGFileSystem().base[0].split('root://eosuser.cern.ch')[-1]}/analysis_httcp/cf.PlotVariables1D/{config.campaign.name}/QCD"
        if not os.path.exists(path):
            os.makedirs(path)
        ratio_path = f"{path}/Ratio"
        if not os.path.exists(ratio_path):
            os.makedirs(ratio_path)
                
        for gidx, group_name in enumerate(complete_groups):
            group = qcd_groups[group_name]
   
            # Get the corresponding histograms of the id, if not present, create a zeroed histogram
            get_hist = lambda h, region_name: (
                h[{"category": hist.loc(group[region_name].id)}]
                if group[region_name].id in h.axes["category"]
                else hist.Hist(*[axis for axis in (h[{"category": [0]}] * 0).axes if axis.name != 'category'])
            ) 

            # Get the corresponding histograms and convert them to number objects,
            os_noniso_mc  = hist_to_num(get_hist(mc_hist, "os_noniso"), "os_noniso_mc")
            os_noniso_data = hist_to_num(get_hist(data_hist, "os_noniso"), "os_noniso_data")

            ## DATA - MC of region C (FF are already apply to them)
            fake_hist = os_noniso_data - os_noniso_mc

            # combine uncertainties and store values in bare arrays
            fake_hist_values = fake_hist()
            fake_hist_variances = fake_hist(sn.UP, sn.ALL, unc=True)**2

            # Guaranty positive values of fake_hist
            neg_int_mask = fake_hist_values <= 0
            fake_hist_values[neg_int_mask] = 1e-5
            fake_hist_variances[neg_int_mask] = 0

            ## Use fake_hist as qcd histogram for category D (os_iso)
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    qcd_hist.view().value[cat_index, ...] = fake_hist_values
                    qcd_hist.view().variance[cat_index, ...] = fake_hist_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {mc_hist} "
                    f"for category {group.os_iso}",
                )

            if type_extrapolation != "CD":
                # Save tne qcd histogram in a pickle file
                hname = qcd_hist.axes[2].name
                #path = "/eos/user/o/oponcet2/analysis/CP_dev/analysis_httcp/cf.PlotVariables1D/QCD"
                # Ensure the folder exists
                #if not os.path.exists(path):
                #    os.makedirs(path)
                with open(f"{path}/qcd_{hname}_{group_name}.pkl", "wb") as f:
                    pickle.dump(qcd_hist, f)

                # create a hist clone of the data_hist
                ratio_hist = data_hist.copy()

                # calultate sum_mc_hist
                mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
                mc_hist_sum = sum(mc_hists[1:], mc_hists[0].copy())

                # For inclusive region
                mc_hist_incl = mc_hist_incl + mc_hist_sum.copy()
                data_hist_incl = data_hist_incl + data_hist.copy()
                
                os_iso_mc  = hist_to_num(get_hist(mc_hist_sum, "os_iso"), "os_iso_mc")
                os_iso_data = hist_to_num(get_hist(data_hist, "os_iso"), "os_iso_data")
                
                # Calucate the DATA/MC ratio
                ratio = os_iso_data/os_iso_mc

                # total MC
                os_iso_mc_incl   = (os_iso_mc + os_iso_mc_incl) if gidx > 0 else os_iso_mc
                os_iso_data_incl = (os_iso_data + os_iso_data_incl) if gidx > 0 else os_iso_data
                
                # combine uncertainties and store values in bare arrays
                ratio_hist_values = ratio()
                ratio_hist_variances = ratio(sn.UP, sn.ALL, unc=True)**2
                
                
                # Guaranty positive values of fake_hist
                neg_int_mask = ratio_hist_values <= 0
                ratio_hist_values[neg_int_mask] = 1e-5
                ratio_hist_variances[neg_int_mask] = 0

                ## Use fake_hist as qcd histogram for category D (os_iso)
                cat_axis = qcd_hist.axes["category"]
                for cat_index in range(cat_axis.size):
                    if cat_axis.value(cat_index) == group.os_iso.id:
                        ratio_hist.view().value[cat_index, ...] = ratio_hist_values
                        ratio_hist.view().variance[cat_index, ...] = ratio_hist_variances
                        break
                    else:
                        raise RuntimeError(
                            f"could not find index of bin on 'category' axis of qcd histogram {mc_hist} "
                            f"for category {group.os_iso}",
                        )
            
                #path = f"{path}/Ratio"
                # save the ratio in a pickle file
                with open(f"{ratio_path}/ratio_{hname}_{group_name}_{group.os_iso.id}.pkl", "wb") as f:
                    pickle.dump(ratio_hist, f)
                    
        if type_extrapolation != "CD":
            # Save the inclusive histograms
            incl_ratio = os_iso_data_incl/os_iso_mc_incl

            incl_ratio_hist_values = incl_ratio()
            incl_ratio_hist_variances = incl_ratio(sn.UP, sn.ALL, unc=True)**2

            incl_ratio_hist = mc_hist_incl.copy().reset()
            incl_ratio_hist.view().value[0, ...] = incl_ratio_hist_values
            incl_ratio_hist.view().variance[0, ...] = incl_ratio_hist_variances
        
            with open(f"{ratio_path}/RATIO_{hname}_inclusive.pkl", "wb") as f:
                pickle.dump(incl_ratio_hist, f) #data_hist_incl, f)
          
        return hists
    
    

    config.x.hist_hooks = {
        "produce_fake_factor": produce_fake_factor,
        "extrapolate_fake": extrapolate_fake
    }
