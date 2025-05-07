import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import comb
import matplotlib.cm as cm


def algo1(noeuds, degre, indice_intervalle):
    S = np.eye(1)
    for k in range(1, degre+1):
        debut = max(0, indice_intervalle - k)
        fin = min(len(noeuds), indice_intervalle + k + 2)
        noeuds_locaux = noeuds[debut:fin]
        
        A = np.zeros((k, k+1))
        B = np.zeros((k, k+1))
        
        for l in range(k):
            idx = l + 1
            denominateur = noeuds_locaux[idx+k] - noeuds_locaux[idx]
            alpha = (noeuds_locaux[k] - noeuds_locaux[idx]) / denominateur if denominateur != 0 else 0.0
            beta = (noeuds_locaux[k+1] - noeuds_locaux[idx]) / denominateur if denominateur != 0 else 0.0
            
            A[l, l] = 1 - alpha
            A[l, l+1] = alpha
            B[l, l] = 1 - beta
            B[l, l+1] = beta
        
        S_haut = S @ A
        S_bas = S[-1:] @ B
        S = np.vstack([S_haut, S_bas])
    
    return S

def bernstein(i, degre, t):
    return comb(degre, i) * (t**i) * (1 - t)**(degre-i)

def evaluer_courbe_bezier(points_controle, poids, degre, nb_points=100):
    t = np.linspace(0, 1, nb_points)
    courbe = np.zeros((nb_points, points_controle.shape[1]))

    for i in range(degre+1):
        B = bernstein(i, degre, t)
        WB = poids[i] * B
        contribution = WB[:, np.newaxis] * points_controle[i]
        courbe += contribution

    denominateur = np.zeros(nb_points)

    for i in range(degre+1):
        B = bernstein(i, degre, t)
        WB = poids[i] * B
        denominateur += WB

    courbe = courbe / denominateur[:, np.newaxis]
    return courbe

def eval_courbe_nurbs(noeuds, points_controle, poids, degre, nb_points=300):
    u_min = noeuds[degre]
    u_max = noeuds[-degre-1]
    u_vals = np.linspace(u_min, u_max, nb_points)
    courbe = np.zeros((nb_points, points_controle.shape[1]))
    
    for idx, u in enumerate(u_vals):
        numerateur = np.zeros(points_controle.shape[1])
        denominateur = 0.0
        for i in range(len(points_controle)):
            N = R_de_Boor(i, degre, noeuds, u)
            numerateur += poids[i] * N * points_controle[i]
            denominateur += poids[i] * N
        courbe[idx] = numerateur / denominateur
    return courbe


def eval_surface_nurbs(noeuds_u, noeuds_v, points_controle, poids, degre_u, degre_v, nb_points_u=50, nb_points_v=50):
    u_min = noeuds[degre_u]
    u_max = noeuds[-degre_u-1]
    u_vals = np.linspace(u_min, u_max, nb_points_u)

    v_min = noeuds[degre_v]
    v_max = noeuds[-degre_v-1]
    v_vals = np.linspace(v_min, v_max, nb_points_v)

    surface = np.zeros((nb_points_u, nb_points_v, 3))

    for iu, u in enumerate(u_vals):
        for iv, v in enumerate(v_vals):
            numerateur = np.zeros(3)
            denominateur = 0.0
            for i in range(points_controle.shape[0]):
                Ni = R_de_Boor_surface(i, degre_u, noeuds_u, u)
                for j in range(points_controle.shape[1]):
                    Mj = R_de_Boor_surface(j, degre_v, noeuds_v, v)
                    poids_ij = poids[i, j]
                    NMi_w = Ni * Mj * poids_ij
                    numerateur += NMi_w * points_controle[i, j]
                    denominateur += NMi_w
            surface[iu, iv] = numerateur / denominateur if denominateur != 0 else numerateur
    return surface


def R_de_Boor(i, degre, noeuds, u):
    n = len(noeuds) - 1
    if degre == 0:
        if i >= n:
            return 0.0
        return 1.0 if noeuds[i] <= u < noeuds[i+1] else 0.0
    terme1 = 0.0
    terme2 = 0.0
    if (i + degre) < n:
        denom1 = noeuds[i+degre] - noeuds[i]
        if denom1 != 0:
            terme1 = (u - noeuds[i]) / denom1 * R_de_Boor(i, degre-1, noeuds, u)
    if (i + degre + 1) < n:
        denom2 = noeuds[i+degre+1] - noeuds[i+1]
        if denom2 != 0:
            terme2 = (noeuds[i+degre+1] - u) / denom2 * R_de_Boor(i+1, degre-1, noeuds, u)
    return terme1 + terme2

#on l'ajoute pour les surfaces specifiquement
def R_de_Boor_surface(i, degre, noeuds, parametre):
    n = len(noeuds) - 1
    if degre == 0:
        if i >= n:
            return 0.0
        return 1.0 if noeuds[i] <= parametre < noeuds[i+1] else 0.0
    terme1 = 0.0
    terme2 = 0.0
    if (i + degre) < n:
        denom1 = noeuds[i+degre] - noeuds[i]
        if denom1 != 0:
            terme1 = (parametre - noeuds[i]) / denom1 * R_de_Boor_surface(i, degre-1, noeuds, parametre)
    if (i + degre + 1) < n:
        denom2 = noeuds[i+degre+1] - noeuds[i+1]
        if denom2 != 0:
            terme2 = (noeuds[i+degre+1] - parametre) / denom2 * R_de_Boor_surface(i+1, degre-1, noeuds, parametre)
    return terme1 + terme2



noeuds = [0, 0, 0, 0, 1/5, 2/5, 2/5, 3/5, 3/5, 3/5, 1, 1, 1, 1]

degre = 3

points_controle = np.array([
    [0, 6, 0],
    [1, 10, 0],
    [5, 12, 0],
    [9, 0, 0],
    [8, 3, 0],
    [5, 1.5, 0],
    [0, 0, 0],
    [2, -2, 0],
    [8, -2, 0],
    [10, 0, 0]
])

poids = np.array([1, 2, 2, 1, 0.5, 0.5, 1, 1, 2, 1])



segments_bezier = []

for i in range(degre, len(noeuds) - degre - 1):
    if noeuds[i] != noeuds[i+1]:
        premier = i - degre
        dernier = i
        if premier < 0 or dernier >= len(points_controle):
            continue
        
        S = algo1(noeuds, degre, i)
        print(f"Matrice S pour l'intervalle {i} :")
        print(S)
        print()
        
        pts_locaux = points_controle[premier:dernier+1]
        poids_locaux = poids[premier:dernier+1]
        
        pointC_pondere = poids_locaux[:, np.newaxis] * pts_locaux
        bezier_PP = S @ pointC_pondere
        bezier_poids = S @ poids_locaux
        
        bezier_pts = bezier_PP / bezier_poids[:, np.newaxis]
        
        courbe = evaluer_courbe_bezier(bezier_pts, bezier_poids, degre)
        segments_bezier.append(courbe)


courbe_nurbs = eval_courbe_nurbs(noeuds, points_controle, poids, degre)



fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(111, projection='3d')
couleurs = cm.get_cmap('tab10', len(segments_bezier))

for idx, segment in enumerate(segments_bezier):
    ax.plot(segment[:,0], segment[:,1], segment[:,2], color=couleurs(idx), label=f'Edge {idx+1}')

# NURBS brute
ax.plot(courbe_nurbs[:,0], courbe_nurbs[:,1], courbe_nurbs[:,2], 'k--', linewidth=2, label='NURBS curve ')

# Points de contrôle
ax.plot(points_controle[:,0], points_controle[:,1], points_controle[:,2], 'ks--', label="Control points")

u_vals = np.linspace(noeuds[degre], noeuds[-degre-1], 1000)
x_vals = []
y_vals = []
z_vals = []
for u in noeuds[degre:-degre]:
    numerateur = np.zeros(points_controle.shape[1])
    denominateur = 0.0
    for i in range(len(points_controle)):
        N = R_de_Boor(i, degre, noeuds, u)
        numerateur += poids[i] * N * points_controle[i]
        denominateur += poids[i] * N
    pt = numerateur / denominateur
    x_vals.append(pt[0])
    y_vals.append(pt[1])
    z_vals.append(pt[2])

ax.scatter(x_vals, y_vals, z_vals, color='red', s=50, label="Bézier points")

ax.set_title("NURBS to Bézier")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.legend()
ax.grid(True)
plt.show()
