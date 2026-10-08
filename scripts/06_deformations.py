"""First-order deformation dimensions (an LU invariant; numerical with a large spectral gap).

dim = nullity of the Jacobian of (rho_A - I/36) over the three cuts  -  rank of the local-unitary orbit tangent (141).
Expected:  psi and the real classes on the gluing locus (1, 2, 9-12): 14 (= 11 face directions + 3 branch directions);
the other real classes (3-8, 13-20): 11;
generic face point: 11;  generic support-72 member and the uniform states: 10;
golden state of Rather et al. (2022): 4;  two-unitary Hadamard H(0) of Bruzda-Zyczkowski (2024): 3.
Different values prove LU inequivalence (even up to party permutations); in particular none of our states is
LU-equivalent to the golden state or to H(0).
"""
import numpy as np

from _common import check, done
from ame46real import core, family as fam
from ame46real.realpoints import Face
from ame46real.jacobian import deformation_dim

rng = np.random.default_rng(7)
psi = core.psi()
states = {"psi": (core.psi(), 14)}
faces = {s: Face(s) for s in "+-"}
for c in core.load_json("real_classes.json")["classes"][1:]:
    F = faces[c["face"]]
    states[f"real class {c['cls']} ({c['face']}, flip {' '.join(c['flip']) or '-'})"] = (F.state(F.flips(c["flip"])),
                                                                                        c["deformation_dim"])
states["face, random theta"] = (fam.psi_face(rng.uniform(0, 2 * np.pi, 11)), 11)
states["support 72, sigma=+1, random"] = (fam.psi72(1, rng.uniform(-1, 1), rng.uniform(0, 2 * np.pi, 9)), 10)
states["support 72, sigma=-1, random"] = (fam.psi72(-1, rng.uniform(-1, 1), rng.uniform(0, 2 * np.pi, 9)), 10)
states["uniform psi_+(0, 0)"] = (fam.psi72(1, 0.0, np.zeros(9)), 10)
golden = np.load(core.DATA + "/reference/golden_rbb22_U36.npy").ravel()
states["golden (Rather et al. 2022)"] = (golden / np.linalg.norm(golden), 4)
B = np.array(core.load_json("reference/bz24_H0_B.json")["B"])
states["H(0) (Bruzda-Zyczkowski 2024)"] = (np.exp(1j * np.pi * B / 3).ravel() / 36, 3)
print(f"{'state':<40} {'AME err':>8} {'dim':>4} {'nullity':>8}  gap (largest zero / smallest nonzero eigenvalue)")
for name, (v, expected) in states.items():
    dim, nul, rT, gap = deformation_dim(v)
    print(f"{name:<40} {core.ame_error(v):8.1e} {dim:>4} {nul:>8}  {gap[0]:.1e} / {gap[1]:.1e}", flush=True)
    check(f"deformation dim of {name} = {expected}", dim == expected and rT == 141 and gap[1] > 1e6 * gap[0])
done()
