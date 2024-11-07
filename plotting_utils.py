# type: ignore
import networkx as nx
import numpy as np


def draw_graph(
    gs, ax, manual_relabels=None, cut=0.3, pos=None, scale=1, label_rotation=15
):
    if manual_relabels is None:
        manual_relabels = {}
    gs_red = gs.edge_subgraph(
        [(u, v) for u, v, d in gs.edges(data=True) if np.abs(d["weight"] / scale) > cut]
    )

    if pos is None:
        pos = nx.shell_layout(gs_red)
    nx.draw_networkx_nodes(
        gs_red,
        pos=pos,
        ax=ax,
        alpha=1.0,
        node_size=[
            50 * np.sqrt(np.abs(gs.in_degree(n, weight="weight") / scale))
            for n in gs_red.nodes
        ],
        node_color=[
            "#aaaaff" if gs.in_degree(n, weight="weight") > 0 else "#ffaaaa"
            for n in gs_red.nodes
        ],
    )

    # draw node labels by hand becasue nx.draw_networkx_labels can't rotate labels
    for label, (x, y) in pos.items():
        if label not in gs_red.nodes():
            continue
        if label in manual_relabels:
            lab = manual_relabels[label]
        else:
            lab = label

        if label_rotation is None:
            theta = np.arctan(y / x) * 180 / np.pi
        else:
            theta = label_rotation
        ax.text(
            x,
            y,
            lab,
            fontsize=12,
            ha="center",
            va="center",
            rotation=theta,
        )

    nx.draw_networkx_edges(
        gs_red,
        pos=pos,
        edge_color=[
            "#aaaaff" if d["weight"] > 0 else "#ffaaaa"
            for _, _, d in gs_red.edges(data=True)
        ],
        alpha=[
            min((np.abs(d["weight"] / scale) / 100) ** (1 / 4) + 0.2, 1)
            for _, _, d in gs_red.edges(data=True)
        ],
        connectionstyle="arc3, rad = 0.1",
        arrowsize=20,
        width=[
            np.log(np.abs(d["weight"] / scale)) + 3
            for _, _, d in gs_red.edges(data=True)
        ],
        ax=ax,
    )
    ax.set_frame_on(False)


def plot_impact_graph(
    ax,
    model,
    fva_healthy,
    fva_phm2,
    amgs,
    manual_relabels,
    cut=0.1,
    sort=True,
    highlights=None,
):
    diffs = fva_phm2 - fva_healthy
    if sort:
        diff_s = sorted(
            diffs.iterrows(),
            key=lambda x: -np.abs(
                (fva_healthy.loc[x[0], "maximum"] + fva_healthy.loc[x[0], "minimum"])
                / 2
                - (fva_phm2.loc[x[0], "maximum"] + fva_phm2.loc[x[0], "minimum"]) / 2
            )
            if x[0] != "PHM2_prodrxn_VN"
            else 1,
        )
    else:
        diff_s = diffs.iterrows()
    barw = 0.85
    ylabs = []
    scales = []
    row = 0
    for i, (dind, (dmin, dmax)) in enumerate(diff_s):
        r_healthy = fva_healthy[fva_healthy.index == dind]
        if r_healthy.empty:
            continue
        hmin = float(r_healthy["minimum"].iloc[0])
        hmax = float(r_healthy["maximum"].iloc[0])
        hwidth = hmax - hmin

        r_infected = fva_phm2[fva_phm2.index == dind]
        imin = float(r_infected["minimum"].iloc[0])
        imax = float(r_infected["maximum"].iloc[0])
        iwidth = imax - imin

        normalization = max([np.abs(hmin), np.abs(hmax), np.abs(imin), np.abs(imax)])

        if normalization < cut:
            continue

        hmin /= normalization
        hmax /= normalization
        hwidth /= normalization
        imin /= normalization
        imax /= normalization
        iwidth /= normalization
        dmin /= normalization
        dmax /= normalization

        hmid = (hmin + hmax) / 2
        imid = (imin + imax) / 2

        if np.abs(hmid - imid) < cut:
            continue

        ax.broken_barh(
            [(hmin, hwidth)],
            (-row - 0.5, barw),
            color="b",
            alpha=0.2,
        )
        ax.plot(hmid, -row, "bo", alpha=0.5)

        ax.broken_barh(
            [(imin, iwidth)],
            (-row - 0.5, barw),
            color="r",
            alpha=0.2,
        )
        ax.plot(imid, -row, "ro", alpha=0.5)

        if highlights and dind in highlights:
            lb, ub = highlights[dind]
            hlb = lb / normalization
            hub = ub / normalization
            ax.broken_barh(
                [(hlb, hub - hlb)],
                (-row - 0.5, barw),
                color="k",
                alpha=0.25,
            )
            ax.plot((hlb + hub) / 2, -row, "ko", alpha=0.5)

        row += 1
        ss = model.reactions.get_by_id(dind).subsystem
        # nm = healthy.reactions.get_by_id(dind).name
        nm = dind
        if not ss:
            ss = "Unspecified"
        if ss in manual_relabels:
            ss = manual_relabels[ss]
        ylabs.append(f"{'*' if nm in amgs else ''}{nm} ({ss})")
        if normalization > 0.01:
            scales.append(f"x{normalization:.2f} " + r"mmol gDW$^{-1}$h$^{-1}$")
        else:
            scales.append(f"x{normalization:.2e} " + r"mmol gDW$^{-1}$h$^{-1}$")
    ax.set_yticks([-x for x in range(len(ylabs))])
    ax.set_yticklabels(ylabs)
    ax2 = ax.twinx()
    ax2.set_yticks([-x for x in range(len(ylabs))])
    ax2.set_yticklabels(scales)
    ax.set_ylim([-len(ylabs), 1])
    ax2.set_ylim([-len(ylabs), 1])
    ax.grid(axis="both")
    ax.set_xlabel("FVA Ranges (normalized)")
