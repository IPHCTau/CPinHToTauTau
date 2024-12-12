import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column,remove_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)



@producer(
    uses={"hcand.decayMode", "hcand.genPartFlav",
          "process_id", "channel_id"},
    produces={"process_id"},
    mc_only=True,
)
def split_dy(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    tau_part_flav = {
        "unknown"   : 0,
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "e->tau"    : 3,
        "mu->tau"   : 4,
        "tau_had"   : 5,
        "jet->tau"  : 6,
    }

    ch_etau_id = self.config_inst.get_channel("etau").id
    ch_mutau_id = self.config_inst.get_channel("mutau").id
    ch_tautau_id = self.config_inst.get_channel("tautau").id
    
    h1_DM = events.hcand.decayMode[:,0]
    h2_DM = events.hcand.decayMode[:,1]

    h1_genflv = events.hcand.genPartFlav[:,0]
    h2_genflv = events.hcand.genPartFlav[:,1]

    h2_match_m = (h2_genflv == tau_part_flav["prompt_mu"]) | (h2_genflv == tau_part_flav["mu->tau"])
    h2_match_e = (h2_genflv == tau_part_flav["prompt_e"]) | (h2_genflv == tau_part_flav["e->tau"])
    h2_match_j = (h2_genflv == tau_part_flav["unknown"]) | (h2_genflv == tau_part_flav["jet->tau"]) 
    h2_genuine = h2_genflv == tau_part_flav["tau_had"]

    h1_match_m = (h1_genflv == tau_part_flav["prompt_mu"]) | (h1_genflv == tau_part_flav["mu->tau"])
    h1_match_e = (h1_genflv == tau_part_flav["prompt_e"]) | (h1_genflv == tau_part_flav["e->tau"])
    h1_match_j = (h1_genflv == tau_part_flav["unknown"]) | (h1_genflv == tau_part_flav["jet->tau"]) 
    h1_genuine = h1_genflv == tau_part_flav["tau_had"]
    
    
    # N.B. For tau-tau channel, lets assume the leading tau is true
    # hardcoded : BAD
    z_j2tau_proc_id = 51097 # if the hcand2 i.e. tauh is from jet or not
    z_t2tau_proc_id = 51098 # if the hcand2 i.e. tauh is true or not
    z_l2tau_proc_id = 51099 # if the hcand2 i.e. tauh is from e/mu or not

    z_t2tau_mask = ak.where(events.channel_id == ch_tautau_id, h1_genuine, h2_genuine)
    z_l2tau_mask = ak.where(events.channel_id == ch_tautau_id, (h1_match_m | h1_match_e), (h2_match_m | h2_match_e))
    z_j2tau_mask = ak.where(events.channel_id == ch_tautau_id, h1_match_j, h2_match_j)
    
    process_id = ak.where(z_l2tau_mask, z_l2tau_proc_id, events.process_id)
    process_id = ak.where(z_j2tau_mask, z_j2tau_proc_id, process_id)
    process_id = ak.where(z_t2tau_mask, z_t2tau_proc_id, process_id)

    events = remove_ak_column(events, "process_id")
    events = set_ak_column(events, "process_id", process_id, value_type=np.int64)

    return events
