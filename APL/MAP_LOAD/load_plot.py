import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, RegularPolygon
from matplotlib.legend_handler import HandlerPatch

flat = 98.5
a = flat / math.sqrt(3.0)

def hex_cells(radius):
    cells = []
    for q in range(-radius, radius + 1):
        rmin = max(-radius, -q - radius)
        rmax = min(radius, -q + radius)
        for r in range(rmin, rmax + 1):
            cells.append((q, r))
    return cells

def hex_ring(radius):
    """Координаты одного кольца в правильном порядке обхода."""
    if radius == 0:
        return [(0, 0)]

    results = []
    q, r = 0, -radius

    directions = [
        (-1, 1),
        (-1, 0),
        (0, 1),
        (1, -1),
        (1, 0),
        (0, -1),
    ]

    for dq, dr in directions:
        for _ in range(radius):
            results.append((q, r))
            q += dq
            r += dr

    return results

cells = hex_cells(6)

corners = {
    (6, 0),
    (0, 6),
    (-6, 6),
    (-6, 0),
    (0, -6),
    (6, -6),
}
cells = [c for c in cells if c not in corners]

print("Количество ТВС:", len(cells))

def axial_to_xy_flat(q, r, side):
    x = 1.5 * side * q
    y = math.sqrt(3) * side * (r + q / 2.0)
    return x, y

def hex_vertices_flat(xc, yc, side):
    angles_deg = [0, 60, 120, 180, 240, 300]
    verts = []
    for ang in angles_deg:
        t = math.radians(ang)
        x = xc + side * math.cos(t)
        y = yc + side * math.sin(t)
        verts.append((x, y))
    return verts

def hex_distance(q, r):
    s = -q - r
    return max(abs(q), abs(r), abs(s))

values = []
for q, r in cells:
    d = hex_distance(q, r)
    if d <= 1:
        val = 1
    elif d <= 3:
        val = 2
    elif d == 4:
        val = 3
    elif d == 5:
        val = 4
    else:
        val = 5
    values.append(val)

color_map = {
    1: "#ffffff",
    2: "#d9d9d9",
    3: "#a6a6a6",
    4: "#737373",
    5: "#404040",
}

enrichment_map = {
    1: 20,
    2: 25,
    3: 30,
    4: 35,
    5: 40,
}

# -----------------------------
# 6 специальных ТВС в пятом кольце
# три пары: (0,1), (10,11), (20,21)
# -----------------------------
ring5 = hex_ring(5)

black_ring5_cells = {
    (-5, 4), (-5, 5),
    (4, 1),  (5, 0),
    (0, -5), (1, -5),
}

fig, ax = plt.subplots(figsize=(10, 11))

for (q, r), val in zip(cells, values):
    x, y = axial_to_xy_flat(q, r, a)
    verts = hex_vertices_flat(x, y, a)

    poly = Polygon(
        verts,
        closed=True,
        facecolor=color_map[val],
        edgecolor="black",
        linewidth=1.0,
        antialiased=False
    )
    ax.add_patch(poly)

    text_color = "white" if val >= 4 else "black"

    ax.text(
        x, y, f"{val}",
        ha="center", va="center",
        fontsize=10,
        color=text_color
    )

    # черный круг в центре специальных ТВС пятого кольца
    if (q, r) in black_ring5_cells :
        circle = plt.Circle(
            (x, y),
            radius=a * 0.22,
            facecolor="black",
            edgecolor="black",
            zorder=5
        )
        ax.add_patch(circle)

ax.set_aspect("equal", adjustable="box")
ax.autoscale_view()
ax.axis("off")

# -----------------------------
# ЛЕГЕНДА
# -----------------------------
def make_legend_hex(legend, orig_handle, xdescent, ydescent, width, height, fontsize):
    center = (width / 2 - xdescent, height / 2 - ydescent)
    patch = RegularPolygon(
        center,
        numVertices=6,
        radius=min(width, height) / 1.6,
        orientation=math.radians(30)
    )
    patch.set_facecolor(orig_handle.get_facecolor())
    patch.set_edgecolor("black")
    patch.set_linewidth(1.2)
    return patch

legend_handles = [
    RegularPolygon(
        (0, 0),
        numVertices=6,
        radius=6,
        orientation=math.radians(30),
        facecolor=color_map[i],
        edgecolor="black",
        label=f"{i} — {enrichment_map[i]}% обогащения"
    )
    for i in [1, 2, 3, 4, 5]
]

ax.legend(
    handles=legend_handles,
    handler_map={RegularPolygon: HandlerPatch(patch_func=make_legend_hex)},
    loc="upper center",
    bbox_to_anchor=(0.5, -0.03),
    ncol=2,
    frameon=False,
    fontsize=12,
    handlelength=2.5,
    handleheight=2.5,
    handletextpad=1.2,
    labelspacing=1.6,
    columnspacing=2.5,
    borderpad=0.8
)

plt.subplots_adjust(bottom=0.14)
plt.tight_layout()
plt.show()