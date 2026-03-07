# ==================================================
# config.py
# File to resume reduction method and E type usable
# ==================================================

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent / "pmns_factory"))

from pmns_factory.pmns_core.operations.reductions.montgomery_reduction import montgomery_reduction
from pmns_factory.pmns_core.operations.reductions.babai_reduction import babai_rounding_limited_reduction
import pmns_factory.pmns_E_type0 as type0



# method available
METHOD_MONTGOMERY = 0
METHOD_BABAI = 1

# configuration
REDUCTION_CONFIG = {
    METHOD_MONTGOMERY: {
        "py_func": montgomery_reduction,
        "c_method_name": "montgomery_reduction",
        "template_file": "montgomery_template.j2"
    },
    METHOD_BABAI: {
        "py_func": babai_rounding_limited_reduction,
        "c_method_name": "babai_reduction",
        "template_file": "babai_template.j2"
    }
}

# type of external reduction polynomial usable
E_TYPE0 = 0

PMNS_CONFIG = {
    E_TYPE0 : type0,
}