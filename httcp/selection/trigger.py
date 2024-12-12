# coding: utf-8

"""
Trigger selection methods.
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, optional_column as opt

from httcp.util import get_objs_p4

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")


@selector(
    uses={
        "run",
        "TrigObj.id", "TrigObj.pt", "TrigObj.eta", "TrigObj.phi", "TrigObj.filterBits",
    },
    exposed=True,
)
def trigger_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    """
    HLT trigger path selection.
    """
    any_fired = False
    any_fired_all_legs_match = False

    trigger_names                    = []
    trigger_types                    = []
    trigger_ids                      = []
    leg_matched_trigobj_idxs_concat  = []
    leg_matched_trigobjs_concat      = []
    leg_min_pt_concat                = []
    leg_max_eta_concat               = []

    fired_and_all_legs_match_concat  = []
    
    # index of TrigObj's to repeatedly convert masks to indices
    index = ak.local_index(events.TrigObj)
    
    for trigger in self.config_inst.x.triggers:
        # skip the trigger if it does not apply to the dataset
        if not trigger.applies_to_dataset(self.dataset_inst):
            continue

        is_single_e  = trigger.has_tag("single_e")
        is_cross_e   = trigger.has_tag("cross_e_tau")
        is_single_mu = trigger.has_tag("single_mu")
        is_cross_mu  = trigger.has_tag("cross_mu_tau")
        is_cross_tau = trigger.has_tag("cross_tau_tau")
        is_cross_tau_jet = trigger.has_tag("cross_tau_tau_jet")

        trigger_type_temp = ak.Array([""])
        if is_single_e:
            trigger_type_temp = ak.Array(["single_e"])
        elif is_cross_e:
            trigger_type_temp = ak.Array(["cross_e_tau"])
        elif is_single_mu:
            trigger_type_temp =	ak.Array(["single_mu"])
        elif is_cross_mu:
            trigger_type_temp = ak.Array(["cross_mu_tau"])
        elif is_cross_tau:
            trigger_type_temp = ak.Array(["cross_tau_tau"])
        elif is_cross_tau_jet:
            trigger_type_temp = ak.Array(["cross_tau_tau_jet"])
            
        trigger_name_temp = ak.Array([trigger.name])
        trigger_name_array, _ = ak.broadcast_arrays(trigger_name_temp, events.event)
        trigger_names.append(trigger_name_array[:,None])

        trigger_type_array, _ = ak.broadcast_arrays(trigger_type_temp, events.event)
        trigger_types.append(trigger_type_array[:,None])
        
        # get bare decisions
        fired = events.HLT[trigger.hlt_field] == 1

        # apply the run-range
        if self.dataset_inst.is_data:
            if trigger.run_range:
                if trigger.run_range[0] is None: 
                    fired = fired & (events.run < trigger.run_range[1])
                elif trigger.run_range[1] is None: 
                    fired = fired & (events.run > trigger.run_range[0])
                else: 
                    fired = fired & (events.run >= trigger.run_range[0]) & (events.run <= trigger.run_range[1])

        any_fired = any_fired | fired

        
        # get trigger objects for fired events per leg
        leg_matched_trigobj_idxs_shallow = []
        leg_matched_trigobj_idxs         = []
        leg_matched_trigobjs             = []
        leg_min_pt                       = []
        leg_max_eta                      = []
        leg_matched_trigobj_idxs_concat_legs  = None
        leg_matched_trigobjs_concat_legs      = None
        all_legs_match = True

        for leg in trigger.legs:
            # start with a True mask
            leg_mask = abs(events.TrigObj.id) >= 0
            # pdg id selection
            if leg.pdg_id is not None:
                leg_mask = leg_mask & (abs(events.TrigObj.id) == leg.pdg_id)
            # pt cut
            if leg.min_pt is not None:
                min_pt = leg.min_pt * ak.ones_like(events.event)
                leg_min_pt.append(min_pt[:,None])
                leg_mask = leg_mask & (events.TrigObj.pt >= leg.min_pt)
            # eta cut
            if leg.max_abseta is not None:
                max_eta = leg.max_abseta * ak.ones_like(events.event)
                leg_max_eta.append(max_eta[:,None])
                leg_mask = leg_mask & (np.abs(events.TrigObj.eta) <= leg.max_abseta)
            # trigger bits match
            if leg.trigger_bits is not None:
                # OR across bits themselves, AND between all decision in the list
                for bits in leg.trigger_bits:
                    # https://github.com/uhh-cms/hh2bbww/blob/master/hbw/selection/trigger.py#L94
                    leg_mask = leg_mask & ((events.TrigObj.filterBits & bits) == bits)

            leg_matched_trigobj_idxs_shallow.append(index[leg_mask])         # O L D
            leg_matched_trigobj_idxs.append(index[leg_mask][:,None])         # N E W
            leg_matched_trigobjs.append(get_objs_p4(events.TrigObj[index[leg_mask]])[:,None])

            # at least one object must match this leg
            all_legs_match = all_legs_match & ak.any(leg_mask, axis=1)

        # final trigger decision
        fired_and_all_legs_match = fired & all_legs_match
        fired_and_all_legs_match_concat.append(fired_and_all_legs_match[:,None])

        # store all intermediate results for subsequent selectors
        #trigger_data.append((trigger, fired_and_all_legs_match, leg_matched_trigobj_idxs_shallow))

        # store the trigger id
        #ids = ak.where(fired_and_all_legs_match, np.float32(trigger.id), np.float32(np.nan))
        #trigger_ids.append(ak.singletons(ak.nan_to_none(ids)))
        ids = trigger.id * ak.ones_like(events.event)
        ids = ids[:,None]
        trigger_ids.append(ids)
        
        # store the trigger obj indices matched to the trigger legs
        leg_matched_trigobj_idxs_concat_legs = ak.concatenate([*leg_matched_trigobj_idxs], axis=1)
        leg_matched_trigobj_idxs_concat.append(leg_matched_trigobj_idxs_concat_legs[:,None])

        # store the corresponding trigger obj p4s matched to the trigger legs
        leg_matched_trigobjs_concat_legs = ak.concatenate([*leg_matched_trigobjs], axis=1)
        leg_matched_trigobjs_concat.append(leg_matched_trigobjs_concat_legs[:,None])
        
        # store the trigger legs pt threshold mentioned in the triggers description [probably redundant]
        leg_min_pt_concat_legs = ak.concatenate([*leg_min_pt], axis=1)
        leg_min_pt_concat.append(leg_min_pt_concat_legs[:,None])

        # store the trigger legs pt threshold mentioned in the triggers description [probably redundant]
        leg_max_eta_concat_legs = ak.concatenate([*leg_max_eta], axis=1)
        leg_max_eta_concat.append(leg_max_eta_concat_legs[:,None])

        
    fired_and_all_legs_match_concat = ak.concatenate([*fired_and_all_legs_match_concat], axis=1)
    trigger_ids              = ak.concatenate(trigger_ids, axis=1)    
    trigger_names            = ak.concatenate(trigger_names, axis=1)
    trigger_types            = ak.concatenate(trigger_types, axis=1)
    leg_min_pt_concat        = ak.concatenate([*leg_min_pt_concat], axis=1)
    leg_max_eta_concat       = ak.concatenate([*leg_max_eta_concat], axis=1)
    leg_matched_trigobj_idxs = ak.concatenate([*leg_matched_trigobj_idxs_concat], axis=1)
    leg_matched_trigobjs     = ak.concatenate([*leg_matched_trigobjs_concat], axis=1)
    

    # applying the main mask: fired_and_all_legs_match_concat
    # comments: considering 3 events
    # names:
    # [
    #  [ 'HLT_Ele27_WPTight_Gsf' ],
    #  [ 'HLT_Ele27_WPTight_Gsf','HLT_Ele32_WPTight_Gsf','HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1','HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1' ],
    #  [ 'HLT_IsoMu27', 'HLT_IsoMu24' ]
    # ]
    trigger_names_filtered = ak.drop_none(ak.mask(trigger_names,      fired_and_all_legs_match_concat))
    # types:
    # [
    #  [ 'single_e' ],
    #  [ 'single_e','single_e','cross_e_tau','cross_tau_tau' ],
    #  [ 'single_mu','single_mu' ]
    # ]
    trigger_types_filtered = ak.drop_none(ak.mask(trigger_types,      fired_and_all_legs_match_concat))
    # ids:
    # [
    #  [ 111000 ],
    #  [ 111000, 112000, 11151, 15153 ],
    #  [ 131000, 132000 ]
    # ]
    trigger_ids_filtered   = ak.drop_none(ak.mask(trigger_ids,        fired_and_all_legs_match_concat))
    # minpt of both legs combined
    # leg_minpt:
    # [
    #  [ [28.0] ],
    #  [ [28.0], [33.0], [25.0, 35.0], [40.0, 40.0] ], # -- because, last two triggers are cross-triggers
    #  [ [25.0], [25.0] ]
    # ]
    leg_minpt_filtered     = ak.drop_none(ak.mask(leg_min_pt_concat,  fired_and_all_legs_match_concat))
    # for simplicity, save the minpt separated for two legs
    # leg1_minpt:
    # [
    #  [ 28.0 ],
    #  [ 28.0, 33.0, 25.0, 40.0 ],
    #  [ 25.0, 25.0 ]
    # ]
    leg1_minpt_filtered    = ak.firsts(leg_minpt_filtered[:,:,0:1], axis=-1) # 
    # leg2_minpt:
    # will contain None because, single triggers do not have 2nd leg
    # so fill all the None with -1.0
    # [
    #  [ -1.0 ],
    #  [ -1.0, -1.0, 35.0, 40.0 ],
    #  [ -1.0, -1.0 ]
    # ]
    leg2_minpt_filtered    = ak.fill_none(ak.firsts(leg_minpt_filtered[:,:,1:2], axis=-1), -1.0)


    leg_maxeta_filtered     = ak.drop_none(ak.mask(leg_max_eta_concat,  fired_and_all_legs_match_concat))
    leg1_maxeta_filtered    = ak.firsts(leg_maxeta_filtered[:,:,0:1], axis=-1) # 
    leg2_maxeta_filtered    = ak.fill_none(ak.firsts(leg_maxeta_filtered[:,:,1:2], axis=-1), -1.0)

        
    # leg mathced trig obj indices:
    # [
    #  [ [[0]] ],
    #  [ [[0]],   [[0]],   [[0], [2, 37, 38]],   [[2, 37, 38], [2, 37, 38]] ],
    #  [ [[0]],   [[0]] ]
    # ]
    leg_matched_trigobj_idxs_filtered = ak.drop_none(ak.mask(leg_matched_trigobj_idxs, fired_and_all_legs_match_concat))
    # for simplicity, save the matched indices per leg
    # leg 1 mathced trig obj indices:
    # [
    #  [ [0] ],
    #  [ [0],   [0],   [0],   [2, 37, 38] ],
    #  [ [0],   [0] ]
    # ]
    leg1_matched_trigobj_idxs_filtered = ak.firsts(leg_matched_trigobj_idxs_filtered[:,:,0:1], axis=2)
    # leg 2 mathced trig obj indices:
    # [
    #  [ None ],
    #  [ None,   None,   [2, 37, 38],   [2, 37, 38] ],
    #  [ None,   None ]
    # ]
    leg2_matched_trigobj_idxs_filtered = ak.firsts(leg_matched_trigobj_idxs_filtered[:,:,1:2], axis=2)

    
    # leg mathced trig objects p4:
    # [
    #  [ [[p4]] ],
    #  [ [[p4]],   [[p4]],   [[p4], [p4, p4, p4]],   [[p4, p4, p4], [p4, p4, p4]] ],
    #  [ [[p4]],   [[p4]] ]
    # ]
    leg_matched_trigobjs_filtered     = ak.drop_none(ak.mask(leg_matched_trigobjs, fired_and_all_legs_match_concat))
    # for simplicity, save the matched objects per leg
    # leg 1 mathced trig objects:
    # [
    #  [ [p4] ],
    #  [ [p4],   [p4],   [p4],   [p4, p4, p4] ],
    #  [ [p4],   [p4] ]
    # ]
    leg1_matched_trigobjs_filtered = ak.firsts(leg_matched_trigobjs_filtered[:,:,0:1], axis=2)
    # leg 2 mathced trig objects:
    # [
    #  [ None ],
    #  [ None,   None,   [p4, p4, p4],   [p4, p4, p4] ],
    #  [ None,   None ]
    # ]
    leg2_matched_trigobjs_filtered = ak.firsts(leg_matched_trigobjs_filtered[:,:,1:2], axis=2)


    
    # store the fired trigger ids and others
    # e.g.
    #   trigger_names             = [ ['HLT_Ele27_WPTight_Gsf'], [], ['HLT_Ele27_WPTight_Gsf', 'HLT_Ele32_WPTight_Gsf'], ['HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1] ]
    #   trigger_types             = [ [      'single_e'       ], [], [        'single_e'     ,        'single_e'      ], [                          'cross_e_tau'                            ] ]
    #   trigger_ids               = [ [        111000         ], [], [          111000       ,          112000        ], [                              11151                                ] ]
    ##  leg_minpt                 = [ [        [28.0]         ], [], [          [28.0]       ,          [33.0]        ], [                          [25.0, 35.0]                             ] ]
    #   leg1_minpt                = [ [         28.0          ], [], [           28.0        ,           33.0         ], [                               25.0                                ] ]
    #   leg2_minpt                = [ [         -1.0          ], [], [           -1.0        ,           -1.0         ], [                               35.0                                ] ]
    ##  leg_matched_trigobj_idxs  = [ [         [[0]]         ], [], [           [[0]]       ,           [[0]]        ], [                         [[0], [72, 73]]                           ] ]
    ##  leg1_matched_trigobj_idxs = [ [          [0]          ], [], [            [0]        ,            [0]         ], [                                [0]                                ] ]
    ##  leg2_matched_trigobj_idxs = [ [         None          ], [], [           None        ,           None         ], [                              [72,73]                              ] ]
    ##  leg_matched_trigobjs      = [ [        [[p4]]         ], [], [          [[p4]]       ,          [[p4]]        ], [                         [[p4], [p4, p4]]                          ] ]
    #   leg1_matched_trigobjs     = [ [         [p4]          ], [], [           [p4]        ,           [p4]         ], [                               [p4]                                ] ]
    #   leg2_matched_trigobjs     = [ [         None          ], [], [           None        ,           None         ], [                             [p4, p4]                              ] ]

    trigger_data = {
        "trigger_names"  : trigger_names_filtered,
        "trigger_types"  : trigger_types_filtered,
        "trigger_ids"    : trigger_ids_filtered,
        "leg1_minpt"     : leg1_minpt_filtered,
        "leg2_minpt"     : leg2_minpt_filtered,
        "leg1_maxeta"     : leg1_maxeta_filtered,
        "leg2_maxeta"     : leg2_maxeta_filtered,
        "leg1_matched_trigobj_idxs" : leg1_matched_trigobj_idxs_filtered,
        "leg2_matched_trigobj_idxs" : leg2_matched_trigobj_idxs_filtered,
        "leg1_matched_trigobjs" : leg1_matched_trigobjs_filtered,
        "leg2_matched_trigobjs" : leg2_matched_trigobjs_filtered,
    }
    
    return events, SelectionResult(
        steps={
            "trigger": any_fired,
        },
        aux=trigger_data,
    )


@trigger_selection.init
def trigger_selection_init(self: Selector) -> None:
    if getattr(self, "dataset_inst", None) is None:
        return

    # full used columns
    self.uses |= {
        opt(trigger.name)
        for trigger in self.config_inst.x.triggers
        if trigger.applies_to_dataset(self.dataset_inst)
    }
