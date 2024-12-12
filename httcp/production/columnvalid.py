# coding: utf-8
"""                                                                                                                                                                                                               
Column production methods related to higher-level features.
"""
import functools

import law
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.columnar_util import set_ak_column
from columnflow.util import maybe_import
from httcp.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")

logger = law.logger.get_logger(__name__)

@producer(
    uses={
        # nano columns
        "PuppiMET.covXX", "PuppiMET.covXY", "PuppiMET.covYY",
        "PuppiMET.significance",
    },
    produces={
        "PuppiMET.covXX", "PuppiMET.covXY", "PuppiMET.covYY",
        "PuppiMET.significance",
    },
)
def make_column_valid(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    events = set_ak_column(events, "PuppiMET.covXX", ak.nan_to_num(events.PuppiMET.covXX, nan=-9999.9))
    events = set_ak_column(events, "PuppiMET.covXY", ak.nan_to_num(events.PuppiMET.covXY, nan=-9999.9))
    events = set_ak_column(events, "PuppiMET.covYY", ak.nan_to_num(events.PuppiMET.covYY, nan=-9999.9))
    events = set_ak_column(events, "PuppiMET.significance", ak.nan_to_num(events.PuppiMET.significance, nan=-9999.9))
    
    return events
