import law

from typing import Optional
from operator import and_
from functools import reduce
from collections import defaultdict, OrderedDict

from columnflow.selection import Selector, SelectionResult, selector

from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, flat_np_view



logger = law.logger.get_logger(__name__)

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")


@selector(
    uses={
        "cahnnel_id",
        "hcand.*", "hcandprod.*",
    },
    mc_only=True,
)
def find_ghost(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    chid = events.channel_id

    hcands = events.hcand
    hcands_etau = hcands[chid == 1]
    h1_etau = ak.flatten(hcands_etau[:,0:1])
    h2_etau = ak.flatten(hcands_etau[:,1:2])
    
    hcands_mutau = hcands[chid == 2]
    h1_mutau = ak.flatten(hcands_mutau[:,0:1])
    h2_mutau = ak.flatten(hcands_mutau[:,1:2])
    
    hcands_tautau = hcands[chid == 4]
    h1_tautau = ak.flatten(hcands_tautau[:,0:1])
    h2_tautau = ak.flatten(hcands_tautau[:,1:2])

    return events

