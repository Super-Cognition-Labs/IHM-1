#!/usr/bin/env python3
"""The anatomy atlas: whole-body views, organ systems, sections and a peel-away dissection.

Everything here is the canonical anatomy in its own frame (metres, +x left, +y up, +z
anterior), unposed: no plant, no contact, no motion -- the figure shows what the body IS,
not what the scaffold did. Reuses scripts/render_body_3d.py's z-buffered rasteriser.

Licences, enforced in `public()` below and stated on every figure's caption:
  BodyParts3D, (c) The Database Center for Life Science, CC BY 4.0     (body-bp3d-*)
  Z-Anatomy (github.com/Z-Anatomy, rev 702db27), CC BY-SA 4.0          (body-za-*)
The canonical anatomy carries nothing else (audited: 2,229 + 1,767 + 3 synthesized skin
layers + 1 graph). Z-Anatomy's non-commercial inner-ear and kidney parts are excluded by
name here as well, conservatively, even where the source is ambiguous: the ossicles and
renal pelvis are dropped from every figure.

Writes PNG (full resolution) and JPG (web) into docs/figures/atlas/ by default.
"""
import argparse, gzip, json
from pathlib import Path
import importlib.util
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rb3", ROOT / "scripts/render_body_3d.py")
RB = importlib.util.module_from_spec(spec); spec.loader.exec_module(RB)

NC_NAME_KEYS = ("kidney", "renal pelvis", "cochlea", "labyrinth", "vestibul", "semicircular",
                "malleus", "incus", "stapes", "ossicle", "tympan")
SYSTEM_RGB = {
    "skin": (0.87, 0.72, 0.62), "muscular": (0.70, 0.25, 0.23), "skeletal": (0.93, 0.90, 0.82),
    "arterial": (0.80, 0.12, 0.13), "venous": (0.22, 0.32, 0.70), "cardiac": (0.68, 0.13, 0.16),
    "nervous": (0.95, 0.80, 0.25), "lymphatic": (0.40, 0.72, 0.38), "respiratory": (0.88, 0.58, 0.62),
    "digestive": (0.86, 0.58, 0.42), "urinary": (0.78, 0.66, 0.30), "reproductive": (0.80, 0.52, 0.66),
    "endocrine": (0.62, 0.45, 0.75), "sensory": (0.55, 0.60, 0.80), "connective": (0.86, 0.83, 0.70),
    "integumentary": (0.87, 0.72, 0.62),
}
CARTILAGE = (0.72, 0.82, 0.86)
BG_TOP, BG_BOTTOM = np.array([0.975, 0.972, 0.965]), np.array([0.90, 0.895, 0.885])
INK, INK_SOFT = (34, 34, 38), (100, 100, 108)


def public(e):
    if not (e["id"].startswith("body-bp3d-") or e["id"].startswith("body-za-")): return False
    if (e.get("reference_geometry") or {}).get("representation") != "triangular_surface": return False
    if e["role"] in ("surface_region",): return False          # skin-region overlays z-fight the skin
    return not any(k in e["name"].lower() for k in NC_NAME_KEYS)


class Atlas:
    def __init__(self):
        self.entities = [e for e in json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"] if public(e)]
        self.cache = {}
        self.skin = self._skin()

    def geometry(self, e):
        if e["id"] not in self.cache:
            g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
            self.cache[e["id"]] = (np.asarray(g["positions"], float).reshape(-1, 3),
                                   np.asarray(g["indices"], np.int64).reshape(-1, 3))
        return self.cache[e["id"]]

    def _skin(self):
        # The exterior surface the contact work audited (a BodyParts3D skin, CC BY 4.0).
        mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
        e = next(x for x in mech["entities"] if x["role"] == "skin")
        g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
        V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
        ext = np.asarray(json.loads((ROOT / "data/research/engineered_skin_territories/materialization.json").read_text())
                         ["contact_eligible_triangle_ids"], np.int64)
        return V, F[ext]

    def layer(self, pred, rgb=None):
        out = []
        for e in self.entities:
            if not pred(e): continue
            V, F = self.geometry(e)
            col = rgb or (CARTILAGE if e["role"] == "cartilage" else SYSTEM_RGB.get(e["system"], (0.7, 0.7, 0.7)))
            out.append((V, F, col))
        return out


def keep(V, F, test):
    """Triangles whose three vertices all pass `test` (a vertex predicate)."""
    ok = test(V)
    return V, F[ok[F].all(1)]


def clip_scene(scene, test):
    return [(V, keep(V, F, test)[1], c) for V, F, c in scene]


def render(scene, direction, W, H, ss=2, fov=20.0, centre=None, extent=None, up=(0, 1, 0)):
    """Returns (rgb float HxWx3 before downsample, coverage mask) at W*ss x H*ss."""
    scene = [(V, F, c) for V, F, c in scene if len(F)]
    Wi, Hi = W * ss, H * ss
    yy = np.linspace(0, 1, Hi)[:, None, None]
    img = (BG_TOP * (1 - yy) + BG_BOTTOM * yy) * np.ones((Hi, Wi, 3))
    mask = np.zeros(Hi * Wi, bool)
    if not scene: return img, mask.reshape(Hi, Wi)
    d = np.asarray(direction, float); d /= np.linalg.norm(d)
    eye = centre + d * (extent / 2 / np.tan(np.radians(fov / 2)) * 1.05)
    cam = RB.Camera(eye, centre, up=up, fov_deg=fov, aspect=W / H)
    V, F, C, N, off = [], [], [], [], 0
    for Vm, Fm, rgb in scene:
        used = np.unique(Fm); remap = -np.ones(len(Vm), np.int64); remap[used] = np.arange(len(used))
        Vu = Vm[used]; Fu = remap[Fm]
        V.append(Vu); F.append(Fu + off); C.append(np.tile(rgb, (len(Vu), 1))); N.append(RB.vertex_normals(Vu, Fu)); off += len(Vu)
    V, F, C, N = np.vstack(V), np.vstack(F), np.vstack(C), np.vstack(N)
    sx, sy, z, _ = cam.project(V, Wi, Hi)
    r = RB.rasterise(sx, sy, z, F, Wi, Hi, backface=False)
    if r is not None:
        pix, tri, bary, _ = r
        n = (N[F[tri]] * bary[..., None]).sum(1); n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
        view = -cam.forward
        inside = (n @ view) < 0
        n[inside] *= -1                               # two-sided: a cut shows its inside
        key = view * 0.55 + cam.up * 0.60 - cam.right * 0.40; key /= np.linalg.norm(key)
        fill = view * 0.7 + cam.right * 0.55 - cam.up * 0.1; fill /= np.linalg.norm(fill)
        col = (C[F[tri]] * bary[..., None]).sum(1)
        col[inside] *= 0.78                           # cut interiors read darker than surfaces
        half = key + view; half /= np.linalg.norm(half)
        rim = np.clip(1 - np.abs(n @ view), 0, 1) ** 3
        light = 0.30 + 0.62 * np.clip(n @ key, 0, None) + 0.16 * np.clip(n @ fill, 0, None)
        spec_ = 0.09 * np.clip(n @ half, 0, None) ** 30
        out = np.clip(col * light[:, None] + spec_[:, None] - 0.10 * rim[:, None], 0, 1)
        img.reshape(-1, 3)[pix] = out
        mask[pix] = True
    return img, mask.reshape(Hi, Wi)


def finish(img, W, H, ss):
    return Image.fromarray((img.reshape(H, ss, W, ss, 3).mean((1, 3)) * 255).astype(np.uint8))


def ghosted(front, ghost, alpha=0.22):
    """`front` drawn over a faint `ghost` (each an (img, mask) render from the same camera).
    Uncovered pixels of `front` are background, so the ghost blends into the background there."""
    fi, fm = front; gi, gm = ghost
    out = np.where(gm[..., None], (1 - alpha) * fi + alpha * gi, fi)
    return np.where(fm[..., None], fi, out)


def font(size):
    for name in ("DejaVuSans.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        try: return ImageFont.truetype(name, size)
        except OSError: pass
    return ImageFont.load_default()


def label(im, title, sub=None):
    d = ImageDraw.Draw(im); s = im.width / 900
    d.text((int(28 * s), int(22 * s)), title, fill=INK, font=font(int(34 * s)))
    if sub: d.text((int(28 * s), int(66 * s)), sub, fill=INK_SOFT, font=font(int(21 * s)))
    return im


def sheet(panels, caption, gap=0):
    W = sum(p.width for p in panels) + gap * (len(panels) - 1); H = max(p.height for p in panels)
    s = H / 1700
    foot = int(64 * s)
    out = Image.new("RGB", (W, H + foot), tuple(int(v * 255) for v in BG_BOTTOM))
    x = 0
    for p in panels: out.paste(p, (x, 0)); x += p.width + gap
    ImageDraw.Draw(out).text((int(28 * s), H + int(18 * s)), caption, fill=INK_SOFT, font=font(int(20 * s)))
    return out


ATTRIBUTION = ("IHM-1 canonical anatomy, unposed. Meshes: BodyParts3D (c) DBCLS, CC BY 4.0; "
               "Z-Anatomy, CC BY-SA 4.0. Rendered by scripts/render_anatomy_atlas.py.")


def save(im, out, name):
    im.save(out / f"{name}.png")
    web = im.copy(); web.thumbnail((2400, 2400))
    web.save(out / f"{name}.jpg", quality=86, optimize=True)
    print(out / f"{name}.png", im.size, flush=True)


SUPERFICIAL = ("pectoralis major", "deltoid", "trapezius", "latissimus", "external oblique", "gluteus maximus",
               "sartorius", "rectus femoris", "vastus lateralis", "vastus medialis", "gastrocnemius", "biceps brachii",
               "triceps", "sternocleidomastoid", "rectus abdominis", "tensor fasciae", "brachioradialis",
               "gracilis", "adductor longus", "tibialis anterior", "frontalis", "temporalis", "masseter",
               "orbicularis", "platysma", "serratus anterior", "infraspinatus", "teres major", "soleus",
               "biceps femoris", "semitendinosus", "semimembranosus", "flexor carpi", "extensor carpi",
               "palmaris", "pronator teres", "extensor digitorum", "flexor digitorum superficialis")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="docs/figures/atlas")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--W", type=int, default=900)
    ap.add_argument("--H", type=int, default=1700)
    ap.add_argument("--ss", type=int, default=3)
    args = ap.parse_args()
    out = ROOT / args.out; out.mkdir(parents=True, exist_ok=True)
    A = Atlas(); W, H, ss = args.W, args.H, args.ss
    Vs, Fs = A.skin
    centre = (Vs.min(0) + Vs.max(0)) / 2; extent = float((Vs.max(0) - Vs.min(0))[1]) * 1.04
    # The canonical frame is centred (feet at y = -0.865 m), so heights are given above the
    # sole and converted here.
    floor = float(Vs.min(0)[1]); Y = lambda h: floor + h
    want = lambda k: not args.only or k in args.only
    skin = [(Vs, Fs, SYSTEM_RGB["skin"])]
    bones = A.layer(lambda e: e["system"] == "skeletal" or e["role"] in ("rigid_bone", "cartilage"))
    muscles = A.layer(lambda e: e["role"] in ("muscle", "tendon"))
    print(f"public entities {len(A.entities)}; skin {len(Fs):,} tris; bones {len(bones)}; muscles {len(muscles)}", flush=True)

    def view(scene, d, w=None, h=None, **kw):
        return render(scene, d, w or W, h or H, ss, centre=kw.pop("centre", centre), extent=kw.pop("extent", extent), **kw)

    FRONT, BACK, SIDE, THREE = (0, 0, 1), (0, 0, -1), (1, 0, 0), (0.72, 0, 0.69)

    if want("skin"):
        ps = [label(finish(view(skin, d)[0], W, H, ss), t, s) for d, t, s in
              ((FRONT, "Anterior", "exterior skin surface"), (SIDE, "Left lateral", ""), (BACK, "Posterior", ""), (THREE, "Three-quarter", ""))]
        save(sheet(ps, f"The skin: {len(Fs):,} triangles of the exterior surface. " + ATTRIBUTION), out, "01_skin")

    if want("skeleton"):
        ps = [label(finish(view(bones, d)[0], W, H, ss), t, s) for d, t, s in
              ((FRONT, "Skeleton", f"{len(bones)} bones and cartilages"), (SIDE, "Left lateral", ""), (BACK, "Posterior", ""))]
        save(sheet(ps, "The skeleton. " + ATTRIBUTION), out, "02_skeleton")

    if want("muscles"):
        deep = [m for m, e in zip(muscles, [e for e in A.entities if e["role"] in ("muscle", "tendon")])
                if not any(k in e["name"].lower() for k in SUPERFICIAL)]
        ps = [label(finish(view(muscles + bones, FRONT)[0], W, H, ss), "Muscles, anterior", f"{len(muscles)} muscles and tendons"),
              label(finish(view(muscles + bones, BACK)[0], W, H, ss), "Muscles, posterior", ""),
              label(finish(view(deep + bones, FRONT)[0], W, H, ss), "Deep layer, anterior", f"{len(SUPERFICIAL)} superficial groups removed"),
              label(finish(view(deep + bones, BACK)[0], W, H, ss), "Deep layer, posterior", "")]
        save(sheet(ps, "The musculature, superficial and deep. " + ATTRIBUTION), out, "03_muscles")

    if want("systems"):
        ghost = view(skin + bones, FRONT)
        systems = [("Cardiovascular", lambda e: e["system"] in ("arterial", "venous", "cardiac")),
                   ("Nervous", lambda e: e["system"] == "nervous"),
                   ("Respiratory", lambda e: e["system"] == "respiratory"),
                   ("Digestive", lambda e: e["system"] == "digestive"),
                   ("Lymphatic", lambda e: e["system"] == "lymphatic"),
                   ("Urogenital and endocrine", lambda e: e["system"] in ("urinary", "reproductive", "endocrine"))]
        ps = []
        for title, pred in systems:
            layer = A.layer(pred)
            ps.append(label(finish(ghosted(view(layer, FRONT), ghost), W, H, ss), title, f"{len(layer)} structures"))
        save(sheet(ps, "Organ systems, each over a faint skin and skeleton. " + ATTRIBUTION), out, "04_systems")

    everything = lambda: skin + A.layer(lambda e: True)
    if want("sections"):
        cx, cy, cz = centre
        scene = everything()
        sag = clip_scene(scene, lambda V: V[:, 0] <= cx)              # left half removed
        cor = clip_scene(scene, lambda V: V[:, 2] <= cz + 0.005)      # front removed at the mid-coronal plane
        ps = [label(finish(view(sag, (1, 0.12, 0.25))[0], W, H, ss), "Midsagittal section", "left half removed, viewed from the left"),
              label(finish(view(cor, (0.15, 0.10, 1))[0], W, H, ss), "Coronal section", "anterior half removed")]
        for level, name in ((1.30, "Thorax"), (1.05, "Abdomen"), (0.85, "Pelvis")):
            tr = clip_scene(scene, lambda V, l=Y(level): V[:, 1] <= l)
            c = np.array([cx, Y(level) - 0.05, cz]); ps.append(label(finish(view(tr, (0, 1.0, 0.18), centre=c, extent=0.62)[0], W, H, ss),
                                                                 f"Transverse, {name.lower()}", f"cut {level:.2f} m above the sole, viewed from above"))
        save(sheet(ps, "Sections through every layer: skin, muscle, organs, vessels, nerves, bone. " + ATTRIBUTION), out, "05_sections")

    if want("peel"):
        layers = [(skin, -0.09), (muscles, 0.0), (A.layer(lambda e: e["system"] not in ("skeletal", "muscular", "integumentary", "connective") and e["role"] not in ("muscle", "tendon", "rigid_bone", "cartilage")), 0.09), (bones, 9.0)]
        scene = []
        for layer, bound in layers:
            scene += clip_scene(layer, lambda V, b=bound: V[:, 0] <= b)
        whole = label(finish(view(scene, FRONT)[0], W, H, ss), "Dissection", "skin | muscle | organs and vessels | bone")
        tc = np.array([centre[0], Y(1.22), centre[2]])
        torso = label(finish(view(scene, (0.25, 0.08, 1), centre=tc, extent=1.05)[0], W, H, ss), "Dissection, trunk", "")
        tscene = []
        for layer, bound in [(skin, Y(1.38)), (muscles, Y(1.18)), (layers[2][0], Y(0.98)), (bones, -9.0)]:
            tscene += clip_scene(layer, lambda V, b=bound: V[:, 1] >= b)
        steps = label(finish(view(tscene, (0.3, 0.05, 1), centre=tc, extent=1.05)[0], W, H, ss), "Dissection by level", "skin above, then muscle, organs, bone")
        save(sheet([whole, torso, steps], "A peel-away dissection: each layer cut back to show the next. " + ATTRIBUTION), out, "06_dissection")

    if want("head"):
        hc = np.array([centre[0], Y(1.60), centre[2]])
        scene = everything()
        sag = clip_scene(scene, lambda V: (V[:, 0] <= centre[0]) | (V[:, 1] < Y(1.40)))
        ps = [label(finish(view(skin, (0.6, 0.05, 1), h=W, centre=hc, extent=0.34)[0], W, W, ss), "Head", ""),
              label(finish(view(A.layer(lambda e: e["role"] in ("muscle", "tendon", "rigid_bone", "cartilage")), (0.6, 0.05, 1), h=W, centre=hc, extent=0.34)[0], W, W, ss), "Muscles of the head and neck", ""),
              label(finish(view(sag, (1, 0.05, 0.1), h=W, centre=hc, extent=0.34)[0], W, W, ss), "Midsagittal, head", "")]
        save(sheet(ps, "The head and neck. " + ATTRIBUTION), out, "07_head")


if __name__ == "__main__":
    main()
