# proline_multiplier_mapper.py
from typing import Literal
import numpy as np

def map_prob_to_multiplier(
    p: float,
    method: Literal["linear", "sigmoid"] = "linear",
    min_mult: float = 0.5,
    max_mult: float = 2.0,
    gamma: float = 1.0,
    sigmoid_steepness: float = 8.0,
    sigmoid_midpoint: float = 0.5,
) -> float:
    """
    Map a logistic regression probability p (0..1) to a multiplier in [min_mult, max_mult].

    - method "linear": multiplier = min + (max-min) * (p ** gamma)
        - gamma < 1 -> emphasize lower probs, gamma > 1 -> emphasize higher probs
    - method "sigmoid": uses a scaled sigmoid centered at midpoint; safe conservative mapping.

    Returns multiplier rounded to 3 decimals.
    """
    p = float(max(0.0, min(1.0, p)))
    if method == "linear":
        # gamma lets you tune sensitivity (1.0 is straight linear)
        score = p ** float(gamma)
    elif method == "sigmoid":
        # logistic mapping: compress around midpoint
        x = (p - sigmoid_midpoint) * sigmoid_steepness
        s = 1.0 / (1.0 + np.exp(-x))  # 0..1
        # normalize s so s(midpoint)=0.5 maps to 0.5
        score = float(s)
    else:
        raise ValueError("method must be 'linear' or 'sigmoid'")

    mult = min_mult + (max_mult - min_mult) * score
    return round(float(mult), 3)

# -----------------------
# Integration helper
# -----------------------
def model_prob_to_multiplier(
    model,  # trained sklearn-like model with predict_proba(X) method
    X_row,  # feature vector (1d) or dataframe row
    method: Literal["linear","sigmoid"] = "linear",
    **map_kwargs
) -> float:
    """
    Given a trained model and one sample, compute p = P(damaging) and map to multiplier.
    """
    # model.predict_proba expects 2D array; ensure shape
    prob = model.predict_proba([X_row])[0][1]
    return map_prob_to_multiplier(prob, method=method, **map_kwargs)

# -----------------------
# Quick examples
# -----------------------
if __name__ == "__main__":
    # demo: show mapping behavior
    for p in [0.0, 0.05, 0.2, 0.5, 0.8, 0.95, 1.0]:
        print(p, "-> linear:", map_prob_to_multiplier(p, "linear", min_mult=0.5, max_mult=2.0, gamma=1.0),
                 "sigmoid:", map_prob_to_multiplier(p, "sigmoid", min_mult=0.5, max_mult=2.0, sigmoid_steepness=10.0))
