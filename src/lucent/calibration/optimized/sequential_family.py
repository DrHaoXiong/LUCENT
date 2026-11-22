"""Adjust complete-family inputs from sequential.two_stage results."""
import numpy as np
from scipy.stats import false_discovery_control


def adjust_family(family_ids, eligible_ids, results):
    """Require one two-stage result per evaluable positive contrast in the family.

    The caller supplies the frozen full family and the corresponding positive
    contrast IDs. This rejects incomplete runs instead of imputing unrun tests.
    """
    family = tuple(family_ids)
    eligible = tuple(eligible_ids)
    required = set(eligible)
    if not family or len(set(family)) != len(family) or len(required) != len(eligible):
        raise ValueError("Missing family or duplicate identities")
    if not required.issubset(family):
        raise ValueError("Positive contrast outside the family")
    inputs = {}
    for row in results:
        i = row["id"]
        if i in inputs or i not in required:
            raise ValueError("Unexpected or duplicate contrast")
        p = float(row["family_input"])
        if not 0 < p <= 1 or (row["deep"] is None and p != 1):
            raise ValueError("Invalid adjustment input")
        if row["deep"] is not None and p != row["deep"].p:
            raise ValueError("Adjustment input must equal the deep sequential p")
        inputs[i] = p
    if set(inputs) != required:
        raise ValueError("Incomplete positive contrast results")
    values = np.array([inputs.get(i, 1.) for i in family])
    return dict(ids=family, p=values, BY=false_discovery_control(values, method="by"),
                BH=false_discovery_control(values, method="bh"))
