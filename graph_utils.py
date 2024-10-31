from __future__ import annotations

from collections.abc import Mapping
from itertools import combinations
from typing import Any, cast

import cobra  # type: ignore
import networkx as nx  # type: ignore
import requests  # type: ignore
from molmass import Formula


class Reaction:
    id: str
    name: str
    subsystem: str
    upper_bound: float
    lower_bound: float
    products: list[cobra.Metabolite]
    reactants: list[cobra.Metabolite]
    metabolites: dict[cobra.Metabolite, float]


Rxns = list[Reaction]


def diff_graph(
    g1: nx.DiGraph[Any],
    g2: nx.DiGraph[Any],
    weight: str = "weight",
) -> nx.DiGraph[Any]:
    gd: nx.DiGraph[Any] = nx.DiGraph()
    gd.add_nodes_from(g1.nodes)  # type: ignore
    gd.add_nodes_from(g2.nodes)  # type: ignore

    gd.add_edges_from(g1.edges)  # type: ignore
    gd.add_edges_from(g2.edges)  # type: ignore

    for u, v, d in gd.edges(data=True):  # type: ignore
        if g1.has_edge(u, v):
            w1 = g1.edges[u, v].get(weight, 0)
        else:
            w1 = 0

        if g2.has_edge(u, v):
            w2 = g2.edges[u, v].get(weight, 0)
        else:
            w2 = 0
        d[weight] = w1 - w2

    return gd


def _subsystems(model: cobra.Model) -> dict[str, list[str]]:
    subsystems: dict[str, list[str]] = {}
    for rxn in cast(Rxns, model.reactions):  # type: ignore
        if rxn.subsystem in subsystems:
            subsystems[rxn.subsystem].append(rxn.id)
        else:
            subsystems[rxn.subsystem] = [rxn.id]
    return subsystems


def _subsystem_net_metabolites(
    model: cobra.Model,
    fluxes: Mapping[str, float],
    eps: float = 1e-6,
) -> dict[str, dict[str, float]]:
    """Create a mapping from subsystems to a metabolite-flux dictionary."""
    subsystems = _subsystems(model)
    ss_net: dict[str, dict[str, float]] = {}
    for ss, rxs in subsystems.items():
        ss_net[ss] = {}
        for rxn_id in rxs:
            rxn = cast(Reaction, model.reactions.get_by_id(rxn_id))
            for m, stoic in rxn.metabolites.items():
                flux = fluxes[rxn_id] * stoic
                if m.id in ss_net[ss]:
                    ss_net[ss][m.id] += flux
                else:
                    ss_net[ss][m.id] = flux
        for m_id in list(ss_net[ss]):
            if abs(ss_net[ss][m_id]) <= eps:
                del ss_net[ss][m_id]
    return ss_net


def _subsystem_net_metabolites_rev(
    subsystem_net_metabolites: Mapping[str, Mapping[str, float]],
) -> dict[str, dict[str, float]]:
    """Create a mapping from metabolites to a subsystems-flux dictionary."""
    ss_net_rev: dict[str, dict[str, float]] = {}
    for ss, met_dict in subsystem_net_metabolites.items():
        for m, v in met_dict.items():
            if m not in ss_net_rev:
                ss_net_rev[m] = {}
            if ss not in ss_net_rev[m]:
                ss_net_rev[m][ss] = v
            else:
                ss_net_rev[m][ss] += v
    return ss_net_rev


def metabolite_mass_dict(model: cobra.Model) -> dict[str, float | None]:
    masses: dict[str, float | None] = {}
    for met in model.metabolites:  # type: ignore
        if (
            "photon" in met.name.lower()  # type: ignore
            or "ferredoxin" in met.name.lower()  # type: ignore
            or "biomass" in met.name.lower()  # type: ignore
        ):  # type: ignore
            masses[met.id] = 0  # type: ignore
            continue

        formula = cast(str, met.formula)  # type: ignore
        if formula:
            if "R" in formula:
                formula = formula.split("R")[0]
            sv_mass = Formula(formula).isotope.mass
        else:
            name = cast(str, met.name)  # type: ignore
            name = name.removesuffix("[c]")
            name = name.removesuffix("[e]")
            name = name.replace("_", " ")
            name = name.removeprefix("[Protein] ")
            name = name.removeprefix("[Enzyme] ")
            query_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/MolecularWeight/txt"
            ret = [line for line in requests.get(query_url).text.splitlines()]
            if ret[0] != "Status: 404":
                sv_mass = float(ret[0])
            else:
                sv_mass = None
        masses[met.id] = sv_mass  # type: ignore
    return masses


def subsystem_flux_graph(
    model: cobra.Model,
    fluxes: Mapping[str, float],
    eps: float = 1e-6,
    eps_flux_conservation: float = 1e-3,
    masses: dict[str, float | None] | None = None,
) -> nx.DiGraph[str]:
    ss_net = _subsystem_net_metabolites(model, fluxes, eps)
    ss_net_rev = _subsystem_net_metabolites_rev(ss_net)
    if masses is None:
        masses = metabolite_mass_dict(model)
    gs: nx.DiGraph[str] = nx.DiGraph()
    gs.add_nodes_from(ss_net.keys())  # type: ignore
    count = 0
    for s1, s2 in combinations(ss_net, 2):
        s1_to_s2 = 0.0
        s2_to_s1 = 0.0
        for m, sv in ss_net_rev.items():
            sv_mass = masses[m]

            if sv_mass is None:
                count += 1
                print(count, m)
                masses[m] = 0
                sv_mass = 0

            if (s1 not in sv) or (s2 not in sv):
                continue

            # We will double check that the metabolite production is balanced,
            # within numerical tolerance
            tot = sum(x for x in sv.values() if x > 0)
            try:
                assert (
                    tot - sum(-x for x in sv.values() if x < 0) < eps_flux_conservation
                )
            except AssertionError as e:
                print(m, sv)
                print(tot, sum(-x for x in sv.values() if x < 0))
                raise e

            if sv[s1] > 0 and sv[s2] < 0:
                s1_to_s2 = -sv[s1] * sv[s2] / tot
            elif sv[s1] < 0 and sv[s2] > 0:
                s2_to_s1 = -sv[s1] * sv[s2] / tot
        if abs(s1_to_s2) > eps:
            gs.add_edge(s1, s2, weight=sv_mass * s1_to_s2)  # type: ignore
        if abs(s2_to_s1) > eps:
            gs.add_edge(s2, s1, weight=sv_mass * s2_to_s1)  # type: ignore

    gs.add_node("External")  # type: ignore
    for u in gs.nodes():
        if u == "External":
            continue
        deficit = gs.out_degree(u, weight="weight") - gs.in_degree(u, weight="weight")  # type: ignore
        if deficit < 0:
            gs.add_edge(u, "External", weight=-deficit)  # type: ignore
        elif deficit > 0:
            gs.add_edge("External", u, weight=deficit)  # type: ignore

    return gs


def create_metabolome(model: cobra.Model) -> nx.DiGraph[str]:
    metabolome: nx.DiGraph[str] = nx.DiGraph()
    for rxn in cast(Rxns, model.reactions):  # type: ignore
        rxn_id = rxn.id
        subsystem = rxn.subsystem
        if not subsystem:
            subsystem = "unspecified subsystem"
        id = rxn.id
        name = rxn.name
        ids = rxn_id

        # will handle by edges
        backwards = rxn.upper_bound <= 0 and rxn.lower_bound < 0
        reversible = rxn.upper_bound > 0 and rxn.lower_bound < 0

        metabolome.add_node(  # type: ignore
            rxn_id,
            label=rxn_id,
            node_type="Reaction",
            subsystem=subsystem,
            ids=ids,
            backwards=backwards,
            reversible=reversible,
            id=id,
            name=name,
        )

        if reversible:
            metabolome.add_node(  # type: ignore
                rxn_id + "_rev",
                label=rxn_id + "_rev",
                node_type="Reaction",
                subsystem=subsystem,
                ids=ids,
                backwards=backwards,
                reversible=reversible,
                id=id,
                name=name,
            )

        for m, stoic in rxn.metabolites.items():
            metabolome.add_node(m.id, label=m.id, name=m.name, node_type="Metabolite")  # type: ignore
            if backwards:
                if m in rxn.products:
                    metabolome.add_edge(m.id, rxn.id, stoic=abs(stoic))  # type: ignore
                if m in rxn.reactants:
                    metabolome.add_edge(rxn.id, m.id, stoic=abs(stoic))  # type: ignore
            else:  # will apply to reversible also
                if m in rxn.products:
                    metabolome.add_edge(rxn.id, m.id, stoic=abs(stoic))  # type: ignore
                if m in rxn.reactants:
                    metabolome.add_edge(m.id, rxn.id, stoic=abs(stoic))  # type: ignore
            if reversible:
                if m in rxn.products:
                    metabolome.add_edge(m.id, rxn_id + "_rev", stoic=abs(stoic))  # type: ignore
                if m in rxn.reactants:
                    metabolome.add_edge(rxn_id + "_rev", m.id, stoic=abs(stoic))  # type: ignore

    return metabolome


def _normalized_flow_projector(metabolome: nx.DiGraph[str], r1: str, r2: str) -> float:
    """based on the method in:
    Beguerisse-Díaz, Mariano, Gabriel Bosque, Diego Oyarzún, Jesús Picó, and
    Mauricio Barahona. “Flux-Dependent Graphs for Metabolic Networks.” Npj
    Systems Biology and Applications 4, no. 1 (August 14, 2018): 1-14.
    https://doi.org/10.1038/s41540-018-0067-y.

    The basic idea is that the weight of r1->r2 is given by the chance that a
    randomly selected metabolite (stoichiometrically weighted) is produced by r1
    and consumed by r2.
    """
    r1_succ = set(metabolome.successors(r1))
    r2_pred = set(metabolome.predecessors(r2))
    dij: float = 0
    for m in r1_succ & r2_pred:
        wk_prod = metabolome.in_degree(m, weight="stoic")
        wk_cons = metabolome.out_degree(m, weight="stoic")
        s_prod = metabolome.edges[(r1, m)]["stoic"]
        s_cons = metabolome.edges[(m, r2)]["stoic"]
        dij += float(s_prod * s_cons / (wk_cons * wk_prod))  # type: ignore
    return dij


def normalized_flow_graph(metabolome: nx.DiGraph[str]) -> nx.DiGraph[str]:
    top_nodes = set(
        u for u, d in metabolome.nodes(data=True) if d["node_type"] == "Metabolite"
    )
    bot_nodes = set(
        u for u, d in metabolome.nodes(data=True) if d["node_type"] == "Reaction"
    )
    assert top_nodes | bot_nodes == set(metabolome.nodes)

    nfg = nx.bipartite.generic_weighted_projected_graph(  # type: ignore
        metabolome, bot_nodes, weight_function=_normalized_flow_projector
    )

    return cast(nx.DiGraph[str], nfg)
