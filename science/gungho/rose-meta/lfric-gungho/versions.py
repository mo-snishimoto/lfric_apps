import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version31_32 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""

class vn32_t46(MacroUpgrade):
    """Upgrade macro for ticket #46 by Shusuke Nishimoto."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t46"

    def upgrade(self, config, meta_config=None):
        self.remove_setting(config, ["namelist:mixing", "method"])
        self.add_setting(config, ["namelist:mixing", "leonard_tke"], ".true.")
        conf_hash = {}
        conf_list = ["field_names", "enforce_min_value", "min_value"]
        # read list type namelists
        for conf in conf_list:
            conf_val = self.get_setting_value(
                config, ["namelist:transport", conf]
            ).split(",")
            conf_hash[conf] = []
            for value in conf_val:
                if "*" in value:
                    num = int(value.split("*")[0])
                    val = value.split("*")[1]
                    for i in range(num):
                        conf_hash[conf].append(val)
                else:
                    conf_hash[conf].append(value)
        # modify namelist value according to condition
        if (
            "'con_tracer'" in conf_hash["field_names"]
            and "'adv_tracer'" in conf_hash["field_names"]
        ):
            # If both con and adv tracer exist, the former is replaced with
            # pos tracer, and the latter is replaced with gen tracer.
            i = conf_hash["field_names"].index("'con_tracer'")
            j = conf_hash["field_names"].index("'adv_tracer'")
            conf_hash["field_names"][i] = "'pos_tracer'"
            conf_hash["field_names"][j] = "'gen_tracer'"
            conf_hash["enforce_min_value"][i] = ".true."
            conf_hash["enforce_min_value"][j] = ".false."
            conf_hash["min_value"][i] = "0.0"
            conf_hash["min_value"][j] = "-99999999.0"
        elif "'con_tracer'" in conf_hash["field_names"]:
            # If only con tracer exists, it is replaced with pos tracer.
            i = conf_hash["field_names"].index("'con_tracer'")
            conf_hash["field_names"][i] = "'pos_tracer'"
            conf_hash["enforce_min_value"][i] = ".true."
            conf_hash["min_value"][i] = "0.0"
        elif "'adv_tracer'" in conf_hash["field_names"]:
            # If only adv tracer exists, it is replaced with pos tracer.
            i = conf_hash["field_names"].index("'adv_tracer'")
            conf_hash["field_names"][i] = "'pos_tracer'"
            conf_hash["enforce_min_value"][i] = ".true."
            conf_hash["min_value"][i] = "0.0"
        # change namelist value
        for conf in conf_list:
            self.change_setting_value(
                config, ["namelist:transport", conf], ",".join(conf_hash[conf])
            )

        return config, self.reports
