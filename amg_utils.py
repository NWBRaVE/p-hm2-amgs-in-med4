# type: ignore
from cobra.flux_analysis import flux_variability_analysis as fva


def amg_limits(model, amg, direction, tol=0.1, optimum_fraction=1):
    if not direction:
        return
    amg_fva = fva(
        model, [model.reactions.get_by_id(amg)], fraction_of_optimum=optimum_fraction
    )

    lb = amg_fva["minimum"][amg]
    ub = amg_fva["maximum"][amg]
    # print("-" * 80, direction, lb, ub)
    nlb = lb
    nub = ub
    if direction == "raw increase" or direction == "absolute increase":
        nlb = ub - (ub - lb) * tol
    elif direction == "raw decrease":
        nub = lb + (ub - lb) * tol
    elif direction == "absolute decrease":
        if ub <= 0:  # same as raw increase now
            nlb = ub - (ub - lb) * tol
        elif lb >= 0:  # same as raw decrease now
            nub = lb + (ub - lb) * tol
        else:  # ub < 0 and lb > 0
            hw = (ub - lb) * tol / 2
            eu = max(0, hw - ub)  # how far past ub
            el = max(0, lb - hw)  # how far past lb
            nlb = -hw - eu
            nub = hw + el
    else:
        raise ValueError(f"Unknown direction {direction}")

    # floating point errors can arise when the bounds are approximately equal or
    # very close to zero
    if abs(nlb) < 1e-10:
        nlb = 0
    if abs(nub) < 1e-10:
        nub = 0
    if nlb > nub:
        nlb, nub = nub, nlb

    return (nlb, nub)
