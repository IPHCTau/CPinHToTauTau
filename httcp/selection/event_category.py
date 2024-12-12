"""
Exemplary selection methods.
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.util import maybe_import

from httcp.util import transverse_mass
from httcp.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    produces={
        "channel_id",
    },
    exposed=False,
)
def get_categories(
        self: Selector,
        events: ak.Array,
        results: SelectionResult,
        etau_pair_indices: ak.Array,
        mutau_pair_indices: ak.Array,
        tautau_pair_indices: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    # get channels from the config
    ch_etau   = self.config_inst.get_channel("etau")
    ch_mutau  = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")

    false_mask       = (abs(events.event) < 0)

    channel_selections = {
        "cat_is_etau"           : [ch_etau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 2) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 0))],
        "cat_is_mutau"          : [ch_mutau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 0) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 2) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 0))],
        "cat_is_tautau"         : [ch_tautau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 0) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 2))],
        "cat_is_etau_mutau"     : [ch_etau.id + ch_mutau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 2) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 2) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 0))],
        "cat_is_etau_tautau"    : [ch_etau.id + ch_tautau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 2) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 2))], 
        "cat_is_mutau_tautau"   : [ch_mutau.id + ch_tautau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 0) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 2) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 2))], 
    }

    selection_steps  = {}
    channel_id       = np.uint8(1) * false_mask
    for key, val in channel_selections.items():
        selection_steps[key] = val[1]
        channel_id = ak.where(val[1], val[0], channel_id)
        
    # some final type conversions
    channel_id = ak.values_astype(channel_id, np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)

    return events, SelectionResult(
        #steps={
        #    "category_et_mt_or_tt": ((channel_id == ch_etau.id) | (channel_id == ch_mutau.id) | (channel_id == ch_tautau.id)),
        #},
        aux=selection_steps)
