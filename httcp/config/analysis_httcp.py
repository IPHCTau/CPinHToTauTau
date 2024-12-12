# coding: utf-8

"""
Configuration of the CPinHToTauTau analysis.
"""
import law
import order as od

# ------------------------ #
# The main analysis object #
# ------------------------ #
analysis_httcp = ana = od.Analysis(
    name="analysis_httcp",
    id=1,
)

# analysis-global versions
# (see cfg.x.versions below for more info)
ana.x.versions = {}

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.bash_sandboxes = ["$CF_BASE/sandboxes/cf.sh"]
default_sandbox = law.Sandbox.new(law.config.get("analysis", "default_columnar_sandbox"))
if default_sandbox.sandbox_type == "bash" and default_sandbox.name not in ana.x.bash_sandboxes:
    ana.x.bash_sandboxes.append(default_sandbox.name)

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.cmssw_sandboxes = [
    "$CF_BASE/sandboxes/cmssw_default.sh",
]

# config groups for conveniently looping over certain configs
# (used in wrapper_factory)
ana.x.config_groups = {}

# ------------- #
# setup configs #
# ------------- #

# ------------------------------------------------------------- #
#                               Run2                            #
# ------------------------------------------------------------- #
"""
# ===>>> 2018
from httcp.config.config_run2 import add_config as add_config_run2_2018
from cmsdb.campaigns.run2_UL2018_nano_cp_tau_v09 import campaign_run2_UL2018_nano_cp_tau_v09
add_config_run2_2018(
    analysis_httcp,
    campaign_run2_UL2018_nano_cp_tau_v09.copy(),
    config_name=campaign_run2_UL2018_nano_cp_tau_v09.name,
    config_id=int(f"{campaign_run2_UL2018_nano_cp_tau_v09.id}{1}")
)
add_config_run2_2018(
    analysis_httcp,
    campaign_run2_UL2018_nano_cp_tau_v09.copy(),
    config_name=f"{campaign_run2_UL2018_nano_cp_tau_v09.name}_limited",
    config_id=int(f"{campaign_run2_UL2018_nano_cp_tau_v09.id}{2}"),
    limit_dataset_files=1
)
"""
# ------------------------------------------------------------- #
#                               Run3                            #
# ------------------------------------------------------------- #

# ===>>> 2022 PreEE
from httcp.config.config_run3 import add_config as add_config_run3_2022_preEE
from cmsdb.campaigns.run3_2022_preEE_nano_cp_tau_v14 import campaign_run3_2022_preEE_nano_cp_tau_v14
add_config_run3_2022_preEE(
    analysis_httcp,
    campaign_run3_2022_preEE_nano_cp_tau_v14.copy(),
    config_name=campaign_run3_2022_preEE_nano_cp_tau_v14.name,
    config_id=int(f"{campaign_run3_2022_preEE_nano_cp_tau_v14.id}{1}")
)
add_config_run3_2022_preEE(
    analysis_httcp,
    campaign_run3_2022_preEE_nano_cp_tau_v14.copy(),
    config_name=f"{campaign_run3_2022_preEE_nano_cp_tau_v14.name}_limited",
    config_id=int(f"{campaign_run3_2022_preEE_nano_cp_tau_v14.id}{2}"),
    limit_dataset_files=1
)
"""
# ===>>> 2022 PostEE
from httcp.config.config_run3 import add_config as add_config_run3_2022_postEE
from cmsdb.campaigns.run3_2022_postEE_nano_cp_tau_v12 import campaign_run3_2022_postEE_nano_cp_tau_v12
add_config_run3_2022_postEE(
    analysis_httcp,
    campaign_run3_2022_postEE_nano_cp_tau_v12.copy(),
    config_name=campaign_run3_2022_postEE_nano_cp_tau_v12.name,
    config_id=int(f"{campaign_run3_2022_postEE_nano_cp_tau_v12.id}{1}")
)
add_config_run3_2022_postEE(
    analysis_httcp,
    campaign_run3_2022_postEE_nano_cp_tau_v12.copy(),
    config_name=f"{campaign_run3_2022_postEE_nano_cp_tau_v12.name}_limited",
    config_id=int(f"{campaign_run3_2022_postEE_nano_cp_tau_v12.id}{2}"),
    limit_dataset_files=1
)

# ===>>> 2023 PreBPix
from httcp.config.config_run3 import add_config as add_config_run3_2023_preBPix
from cmsdb.campaigns.run3_2023_preBPix_nano_cp_tau_v12 import campaign_run3_2023_preBPix_nano_cp_tau_v12
add_config_run3_2023_preBPix(
    analysis_httcp,
    campaign_run3_2023_preBPix_nano_cp_tau_v12.copy(),
    config_name=campaign_run3_2023_preBPix_nano_cp_tau_v12.name,
    config_id=int(f"{campaign_run3_2023_preBPix_nano_cp_tau_v12.id}{1}")
)
add_config_run3_2023_preBPix(
    analysis_httcp,
    campaign_run3_2023_preBPix_nano_cp_tau_v12.copy(),
    config_name=f"{campaign_run3_2023_preBPix_nano_cp_tau_v12.name}_limited",
    config_id=int(f"{campaign_run3_2023_preBPix_nano_cp_tau_v12.id}{2}"),
    limit_dataset_files=1
)

# ===>>> 2023 PostBPix
from httcp.config.config_run3 import add_config as add_config_run3_2023_postBPix
from cmsdb.campaigns.run3_2023_postBPix_nano_cp_tau_v12 import campaign_run3_2023_postBPix_nano_cp_tau_v12
add_config_run3_2023_postBPix(
    analysis_httcp,
    campaign_run3_2023_postBPix_nano_cp_tau_v12.copy(),
    config_name="campaign_run3_2023_postBPix_nano_cp_tau_v12.name",
    config_id=int(f"{campaign_run3_2023_postBPix_nano_cp_tau_v12.id}{1}")
)
add_config_run3_2023_postBPix(
    analysis_httcp,
    campaign_run3_2023_postBPix_nano_cp_tau_v12.copy(),
    config_name=f"{campaign_run3_2023_postBPix_nano_cp_tau_v12.name}_limited",
    config_id=int(f"{campaign_run3_2023_postBPix_nano_cp_tau_v12.id}{2}"),
    limit_dataset_files=1
)

# ------------------------------------------------------------- #
"""
