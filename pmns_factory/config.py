from pmns_core.operations.reductions import montgomery_reduction
from pmns_core.operations.reductions import babai_rounding_limited_reduction

# method available
METHOD_MONTGOMERY = 0
METHOD_BABAI = 1

# Mapping from names to usable functions
REDUCTION_METHODS = {
    METHOD_MONTGOMERY: montgomery_reduction,
    METHOD_BABAI: babai_rounding_limited_reduction,
}

# Mapping from type of external polynomial usable
E_TYPE_BINOMIAL = 0
E_TYPE_TRINOMIAL = 1

# String of methods name
METHOD_NAMES = {
    METHOD_MONTGOMERY: "Montgomery",
    METHOD_BABAI: "Babai"
}