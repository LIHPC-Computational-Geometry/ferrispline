import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
import matplotlib.cm as cm


def algo1(nodes, degree, indice_intervalle):
    S = np.eye(1)
    for k in range(1, degree + 1):
        debut = max(0, indice_intervalle - k)
        fin = min(len(nodes), indice_intervalle + k + 2)
        nodes_locaux = nodes[debut:fin]

        A = np.zeros((k, k + 1))
        B = np.zeros((k, k + 1))

        for l in range(k):
            idx = l + 1
            denominateur = nodes_locaux[idx + k] - nodes_locaux[idx]
            alpha = (
                (nodes_locaux[k] - nodes_locaux[idx]) / denominateur
                if denominateur != 0
                else 0.0
            )
            beta = (
                (nodes_locaux[k + 1] - nodes_locaux[idx]) / denominateur
                if denominateur != 0
                else 0.0
            )

            A[l, l] = 1 - alpha
            A[l, l + 1] = alpha
            B[l, l] = 1 - beta
            B[l, l + 1] = beta

        S_haut = S @ A
        S_bas = S[-1:] @ B
        S = np.vstack([S_haut, S_bas])

    return S


def bernstein(i, degree, t):
    return comb(degree, i) * (t**i) * (1 - t) ** (degree - i)


def evaluer_courbe_bezier(points_controle, poids, degree, nb_points=100):
    t = np.linspace(0, 1, nb_points)
    courbe = np.zeros((nb_points, points_controle.shape[1]))

    for i in range(degree + 1):
        B = bernstein(i, degree, t)
        WB = poids[i] * B
        contribution = WB[:, np.newaxis] * points_controle[i]
        courbe += contribution

    denominateur = np.zeros(nb_points)

    for i in range(degree + 1):
        B = bernstein(i, degree, t)
        WB = poids[i] * B
        denominateur += WB

    courbe = courbe / denominateur[:, np.newaxis]
    return courbe


def eval_courbe_nurbs(nodes, points_controle, poids, degree, nb_points=300):
    u_min = nodes[degree]
    u_max = nodes[-degree - 1]
    u_vals = np.linspace(u_min, u_max, nb_points)
    courbe = np.zeros((nb_points, points_controle.shape[1]))

    for idx, u in enumerate(u_vals):
        numerateur = np.zeros(points_controle.shape[1])
        denominateur = 0.0
        for i in range(len(points_controle)):
            N = R_de_Boor(i, degree, nodes, u)
            numerateur += poids[i] * N * points_controle[i]
            denominateur += poids[i] * N
        courbe[idx] = numerateur / denominateur
    return courbe


def eval_surface_nurbs(
    nodes_u,
    nodes_v,
    points_controle,
    poids,
    degree_u,
    degree_v,
    nb_points_u=50,
    nb_points_v=50,
):
    u_min = nodes[degree_u]
    u_max = nodes[-degree_u - 1]
    u_vals = np.linspace(u_min, u_max, nb_points_u)

    v_min = nodes[degree_v]
    v_max = nodes[-degree_v - 1]
    v_vals = np.linspace(v_min, v_max, nb_points_v)

    surface = np.zeros((nb_points_u, nb_points_v, 3))

    for iu, u in enumerate(u_vals):
        for iv, v in enumerate(v_vals):
            numerateur = np.zeros(3)
            denominateur = 0.0
            for i in range(points_controle.shape[0]):
                Ni = R_de_Boor_surface(i, degree_u, nodes_u, u)
                for j in range(points_controle.shape[1]):
                    Mj = R_de_Boor_surface(j, degree_v, nodes_v, v)
                    poids_ij = poids[i, j]
                    NMi_w = Ni * Mj * poids_ij
                    numerateur += NMi_w * points_controle[i, j]
                    denominateur += NMi_w
            surface[iu, iv] = (
                numerateur / denominateur if denominateur != 0 else numerateur
            )
    return surface


def R_de_Boor(i, degree, nodes, u):
    n = len(nodes) - 1
    if degree == 0:
        if i >= n:
            return 0.0
        return 1.0 if nodes[i] <= u < nodes[i + 1] else 0.0
    terme1 = 0.0
    terme2 = 0.0
    if (i + degree) < n:
        denom1 = nodes[i + degree] - nodes[i]
        if denom1 != 0:
            terme1 = (u - nodes[i]) / denom1 * R_de_Boor(i, degree - 1, nodes, u)
    if (i + degree + 1) < n:
        denom2 = nodes[i + degree + 1] - nodes[i + 1]
        if denom2 != 0:
            terme2 = (
                (nodes[i + degree + 1] - u)
                / denom2
                * R_de_Boor(i + 1, degree - 1, nodes, u)
            )
    return terme1 + terme2


# on l'ajoute pour les surfaces specifiquement
def R_de_Boor_surface(i, degree, nodes, parametre):
    n = len(nodes) - 1
    if degree == 0:
        if i >= n:
            return 0.0
        return 1.0 if nodes[i] <= parametre < nodes[i + 1] else 0.0
    terme1 = 0.0
    terme2 = 0.0
    if (i + degree) < n:
        denom1 = nodes[i + degree] - nodes[i]
        if denom1 != 0:
            terme1 = (
                (parametre - nodes[i])
                / denom1
                * R_de_Boor_surface(i, degree - 1, nodes, parametre)
            )
    if (i + degree + 1) < n:
        denom2 = nodes[i + degree + 1] - nodes[i + 1]
        if denom2 != 0:
            terme2 = (
                (nodes[i + degree + 1] - parametre)
                / denom2
                * R_de_Boor_surface(i + 1, degree - 1, nodes, parametre)
            )
    return terme1 + terme2


nodes = [0, 0, 0, 0, 1 / 5, 2 / 5, 2 / 5, 3 / 5, 3 / 5, 3 / 5, 1, 1, 1, 1]

degree = 3

points_controle = np.array(
    [
        [0, 6, 0],
        [1, 10, 0],
        [5, 12, 0],
        [9, 0, 0],
        [8, 3, 0],
        [5, 1.5, 0],
        [0, 0, 0],
        [2, -2, 0],
        [8, -2, 0],
        [10, 0, 0],
    ]
)

poids = np.array([1, 2, 2, 1, 0.5, 0.5, 1, 1, 2, 1])


segments_bezier = []

for i in range(degree, len(nodes) - degree - 1):
    if nodes[i] != nodes[i + 1]:
        premier = i - degree
        dernier = i
        if premier < 0 or dernier >= len(points_controle):
            continue

        S = algo1(nodes, degree, i)
        print(f"Matrice S pour l'intervalle {i} :")
        print(S)
        print()

        pts_locaux = points_controle[premier : dernier + 1]
        poids_locaux = poids[premier : dernier + 1]

        pointC_pondere = poids_locaux[:, np.newaxis] * pts_locaux
        bezier_PP = S @ pointC_pondere
        bezier_poids = S @ poids_locaux

        bezier_pts = bezier_PP / bezier_poids[:, np.newaxis]

        courbe = evaluer_courbe_bezier(bezier_pts, bezier_poids, degree)
        segments_bezier.append(courbe)


courbe_nurbs = eval_courbe_nurbs(nodes, points_controle, poids, degree)


fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection="3d")
couleurs = cm.get_cmap("tab10", len(segments_bezier))

for idx, segment in enumerate(segments_bezier):
    ax.plot(
        segment[:, 0],
        segment[:, 1],
        segment[:, 2],
        color=couleurs(idx),
        label=f"Edge {idx+1}",
    )

# NURBS brute
ax.plot(
    courbe_nurbs[:, 0],
    courbe_nurbs[:, 1],
    courbe_nurbs[:, 2],
    "k--",
    linewidth=2,
    label="NURBS curve ",
)

# Points de contrôle
ax.plot(
    points_controle[:, 0],
    points_controle[:, 1],
    points_controle[:, 2],
    "ks--",
    label="Control points",
)

u_vals = np.linspace(nodes[degree], nodes[-degree - 1], 1000)
x_vals = []
y_vals = []
z_vals = []
for u in nodes[degree:-degree]:
    numerateur = np.zeros(points_controle.shape[1])
    denominateur = 0.0
    for i in range(len(points_controle)):
        N = R_de_Boor(i, degree, nodes, u)
        numerateur += poids[i] * N * points_controle[i]
        denominateur += poids[i] * N
    pt = numerateur / denominateur
    x_vals.append(pt[0])
    y_vals.append(pt[1])
    z_vals.append(pt[2])

ax.scatter(x_vals, y_vals, z_vals, color="red", s=50, label="Bézier points")

ax.set_title("NURBS to Bézier")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.legend()
ax.grid(True)
plt.show()
