# type: ignore
from itertools import chain

from cobra.flux_analysis import flux_variability_analysis as fva
from optlang.symbolics import Zero

BOUND_FACTOR = 5
N_LIMIT = -0.86 * BOUND_FACTOR
P_LIMIT = -0.18 * BOUND_FACTOR
C_LIMIT = -5.4 * BOUND_FACTOR
L_LIMIT = -73 * BOUND_FACTOR


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


def amg_impact_by_parsimony(
    model,
    amgs,
    optimum_tolerance=0.5,
    new_amg_penalty=0,
    n_limit=N_LIMIT,
    p_limit=P_LIMIT,
    c_limit=C_LIMIT,
    l_limit=L_LIMIT,
):
    with model as imodel:
        biomass = imodel.slim_optimize()
        imodel.reactions.BIOMASS.lower_bound = biomass * (1 - optimum_tolerance)

        # Do not allow the cell to uptake too much to feed prioritized reactions
        imodel.reactions.AmmoniaEX.bounds = (
            n_limit,
            0,
        )
        imodel.reactions.FAKEOrthophosphateEX.bounds = (
            p_limit,
            0,
        )
        imodel.reactions.HCO3EXcar.bounds = (
            c_limit,
            0,
        )
        imodel.reactions.LightEX.bounds = (
            l_limit,
            0,
        )

        # parsimony only applies to gene-regulated reactions here
        rxn_list = [rxn for rxn in imodel.reactions if rxn.genes or rxn.id in amgs]

        reaction_variables = (
            (rxn.forward_variable, rxn.reverse_variable) for rxn in rxn_list
        )
        reaction_is_doubled = ((rxn.id, rxn.id) for rxn in rxn_list)
        id_var = zip(chain(*reaction_is_doubled), chain(*reaction_variables))

        # Infected
        imodel.objective = imodel.problem.Objective(
            Zero, direction="min", sloppy=True, name="_pfba_objective"
        )

        coeff = {}
        forward = True
        for rid, v in id_var:
            if rid in amgs and amgs[rid] == "raw increase":
                if forward:
                    coeff[v] = 1.0 - abs(new_amg_penalty)
                else:
                    coeff[v] = 1.0 + abs(new_amg_penalty)
            elif rid in amgs and amgs[rid] == "raw decrease":
                if forward:
                    coeff[v] = 1.0 + abs(new_amg_penalty)
                else:
                    coeff[v] = 1.0 - abs(new_amg_penalty)
            elif rid in amgs and amgs[rid] == "absolute increase":
                if forward:
                    coeff[v] = 1.0 - abs(new_amg_penalty)
                else:
                    coeff[v] = 1.0 - abs(new_amg_penalty)
            elif rid in amgs and amgs[rid] == "absolute decrease":
                if forward:
                    coeff[v] = 1.0 + abs(new_amg_penalty)
                else:
                    coeff[v] = 1.0 + abs(new_amg_penalty)
            else:
                coeff[v] = 1.0
            forward = not forward
        imodel.objective.set_linear_coefficients(coeff)
        infected = imodel.optimize()

        return infected
