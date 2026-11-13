# Builds Fig1_composite.html from ../code/Fig1_data.json  (all numbers come from the R simulation)
import json, math, os
here = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(here, "../code/Fig1_data.json")))

TEAL, TEAL_L, ORANGE, SLATE, MUTE, INK, GRAY = "#0F766E", "#99F6E4", "#C2410C", "#334155", "#64748B", "#0F172A", "#94A3B8"
L = d["lines_show"]; S = d["snp_show"]; X = d["X_show"]; G = d["G_show"]; Y = d["y_show"]
B = d["beta_show"]; BB = d["beta_bb_show"]; GV = d["gebv_show"]; GR = d["gebv_rr_show"]
_has1 = [any(r[j] == 1 for r in X) for j in range(len(S))]
_keep = sorted(sorted(range(len(S)), key=lambda j: (not _has1[j], -abs(B[j])))[:8])
S = [S[j] for j in _keep]; X = [[r[j] for j in _keep] for r in X]; B = [B[j] for j in _keep]; BB = [BB[j] for j in _keep]
f = lambda v, k=2: (f"{v:.{k}f}").replace("-", "−")

# ---------- small helpers
def mat(rows, cls="", fmt=lambda v: v, color=None, rowlab=None, collab=None, w=None):
    h = f'<table class="mx {cls}">'
    if collab:
        h += "<tr>" + ("<th></th>" if rowlab else "") + "".join(f'<th class="cl">{c}</th>' for c in collab) + "</tr>"
    for i, r in enumerate(rows):
        h += "<tr>" + (f'<th class="rl">{rowlab[i]}</th>' if rowlab else "")
        for v in (r if isinstance(r, list) else [r]):
            st = f' style="background:{color(v)}"' if color else ""
            h += f"<td{st}>{fmt(v)}</td>"
        h += "</tr>"
    return h + "</table>"

def xcol(v):  return {0: "#F1F5F9", 1: "#99F6E4", 2: "#0F766E"}[v]
def gcol(v):
    if v > 1.5: return "#CCFBF1"
    a = max(-1, min(1, v / 0.45))
    return f"rgba(15,118,110,{a*0.55:.2f})" if a > 0 else f"rgba(194,65,12,{-a*0.45:.2f})"

def rowvec(vals, k=2):
    return '<table class="mx vec"><tr>' + "".join(f"<td>{f(v,k)}</td>" for v in vals) + '<td class="dots" style="width:10px">⋯</td></tr></table>'
def vec(vals, k=2, head=None):
    return '<table class="mx vec">' + "".join(f"<tr><td>{f(v,k)}</td></tr>" for v in vals) + '<tr><td class="dots">⋮</td></tr></table>'

# ---------- Step 1: crosses -> SSD funnel
pairs = d["pairs"]
s1 = '<svg viewBox="0 0 300 158" width="100%" font-family="Inter">'
for k, (a, b) in enumerate(pairs):
    x = 6 + k * 36.5
    s1 += f'<rect x="{x}" y="1" width="33" height="15" rx="4" fill="#fff" stroke="{TEAL_L}"/>'
    s1 += f'<text x="{x+16.5}" y="11.6" text-anchor="middle" font-size="8.5" font-weight="700" fill="{INK}">P{a}×P{b}</text>'
gens = [("F1", 100), ("F2", 50), ("F3", 25), ("F4", 12.5), ("F5", 6.25)]
top_y, h = 28, 19
for i, (g, het) in enumerate(gens):
    y0 = top_y + i * (h + 3); w0 = 236 - i * 22; x0 = 4 + (236 - w0) / 2
    fill = ["#CCFBF1", "#A7F3D0", "#99F6E4", "#5EEAD4", "#2DD4BF"][i]
    s1 += f'<path d="M{x0} {y0} H{x0+w0} L{x0+w0-11} {y0+h} H{x0+11} Z" fill="{fill}"/>'
    s1 += f'<text x="{x0+15}" y="{y0+13}" font-size="10.3" font-weight="800" fill="{INK}">{g}</text>'
    for j in range(8):
        cx = 122 - 3.5 * 11 + j * 11; cy = y0 + 9.5
        s1 += f'<circle cx="{cx}" cy="{cy}" r="3.3" fill="#fff" stroke="{TEAL}" stroke-width="0.9"/>'
        s1 += f'<circle cx="{cx}" cy="{cy}" r="3.3" fill="{TEAL}" opacity="{0.12 + 0.88*(1-het/100):.2f}"/>'
    s1 += f'<text x="296" y="{y0+13}" text-anchor="end" font-size="9.2" font-weight="700" fill="{SLATE}">{het:g}%</text>'
s1 += f'<text x="296" y="{top_y-1}" text-anchor="end" font-size="7.6" fill="{MUTE}">heterozygous</text>'
y_end = top_y + 5 * (h + 3)
s1 += f'<rect x="40" y="{y_end+1}" width="164" height="19" rx="9.5" fill="{TEAL}"/>'
s1 += f'<text x="122" y="{y_end+13.8}" text-anchor="middle" font-size="10.3" font-weight="800" fill="#fff">F5 population: {d["n"]} lines</text>'
s1 += "</svg>"

# ---------- Step 2: marker table
s2 = mat(X, "xm", color=xcol, rowlab=L, collab=[c.replace("SNP", "SNP ") for c in S])
s2 += f'<div class="note">⋯ {d["n"]} lines × {d["n_snp"]} SNPs &nbsp;·&nbsp; 0 = aa, 1 = Aa, 2 = AA</div>'

# ---------- Step 3: PCA
pca = d["pca"]; xs = [p["pc1"] for p in pca]; ys = [p["pc2"] for p in pca]
mnx, mxx, mny, mxy = min(xs), max(xs), min(ys), max(ys)
W3, H3, pad = 230, 118, 26
sx = lambda v: pad + (v - mnx) / (mxx - mnx) * (W3 - pad - 6)
sy = lambda v: H3 - 20 - (v - mny) / (mxy - mny) * (H3 - 58)
s3p = f'<svg viewBox="0 0 {W3} {H3}" width="100%" font-family="Inter">'
s3p += f'<rect x="{pad}" y="6" width="{W3-pad-6}" height="{H3-26}" fill="#F8FAFC" stroke="#E2E8F0"/>'
for p in sorted(pca, key=lambda p: p["train"]):
    if p["train"]:
        s3p += f'<circle cx="{sx(p["pc1"]):.1f}" cy="{sy(p["pc2"]):.1f}" r="2.6" fill="{ORANGE}" stroke="#fff" stroke-width="0.6"/>'
    else:
        s3p += f'<circle cx="{sx(p["pc1"]):.1f}" cy="{sy(p["pc2"]):.1f}" r="2.2" fill="{GRAY}" opacity="0.75"/>'
s3p += f'<text x="{(W3+pad)/2}" y="{H3-6}" text-anchor="middle" font-size="9.2" fill="{SLATE}">PC1 ({d["pv"][0]}%)</text>'
s3p += f'<text x="9" y="{H3/2-8}" text-anchor="middle" font-size="9.2" fill="{SLATE}" transform="rotate(-90 9 {H3/2-8})">PC2 ({d["pv"][1]}%)</text>'
s3p += f'<rect x="{pad+3}" y="7" width="118" height="25" rx="3" fill="#fff" opacity="0.92"/>'
s3p += f'<circle cx="{pad+8}" cy="14" r="2.6" fill="{ORANGE}"/><text x="{pad+14}" y="17" font-size="8.7" font-weight="700" fill="{ORANGE}">training set ({d["n_train"]}), CDmean</text>'
s3p += f'<circle cx="{pad+8}" cy="25" r="2.2" fill="{GRAY}"/><text x="{pad+14}" y="28" font-size="8.7" fill="{SLATE}">untested ({d["n"]-d["n_train"]})</text>'
s3p += "</svg>"
s3g = mat(G, "gm", fmt=lambda v: f(v), color=gcol, rowlab=L, collab=L)

# ---------- Step 4: field
fld = d["field"]
tab = {}
for r in fld: tab.setdefault(r["line"], {})[(r["loc"], r["rep"])] = r["yield"]
cols = [(f"Location {j}", k) for j in (1, 2, 3) for k in (1, 2)]
s4t = '<table class="mx met"><tr><th></th>' + "".join(f'<th class="cl">L{j}<br>R{k}</th>' for j in (1, 2, 3) for k in (1, 2)) + '<th class="cl ym"><i>y</i></th></tr>'
for i, l in enumerate(L):
    s4t += f'<tr><th class="rl">{l}</th>' + "".join(f"<td>{tab[l][c]:.2f}</td>" for c in cols) + f'<td class="ym">{Y[i]:.2f}</td></tr>'
s4t += "</table>"
def plotgrid(x, y, lab):
    g = f'<g transform="translate({x},{y})"><rect width="64" height="44" rx="5" fill="#ECFDF5" stroke="#A7F3D0"/>'
    for r in range(4):
        for c in range(6):
            sh = ["#86EFAC", "#4ADE80", "#22C55E", "#BBF7D0"][(r * 7 + c * 3 + x) % 4]
            g += f'<rect x="{5+c*9.3}" y="{5+r*8.8}" width="7.6" height="7" rx="1.2" fill="{sh}"/>'
    return g + f'<text x="32" y="56" text-anchor="middle" font-size="9.7" font-weight="700" fill="{INK}">{lab}</text></g>'
s4f = '<svg viewBox="0 0 230 66" width="100%" font-family="Inter">' + plotgrid(4, 2, "Location 1") + plotgrid(83, 2, "Location 2") + plotgrid(162, 2, "Location 3") + "</svg>"

# ---------- Step 5 charts: all marker effects
ba = d["beta_all"]; nb = len(ba); mx = max(abs(b["beta"]) for b in ba)
W5, H5 = 300, 50
s5b = f'<svg viewBox="0 0 {W5} {H5}" width="100%" font-family="Inter">'
for k in range(1, 7):
    idx = [i for i, b in enumerate(ba) if b["chr"] == k]
    if k % 2 == 0: s5b += f'<rect x="{6+idx[0]/nb*(W5-10):.1f}" y="4" width="{len(idx)/nb*(W5-10):.1f}" height="{H5-18}" fill="#F1F5F9"/>'
    s5b += f'<text x="{6+(idx[0]+len(idx)/2)/nb*(W5-10):.1f}" y="{H5-3}" text-anchor="middle" font-size="8.0" fill="{MUTE}">Chr {k}</text>'
mid = 4 + (H5 - 18) / 2
s5b += f'<line x1="6" x2="{W5-4}" y1="{mid}" y2="{mid}" stroke="{GRAY}" stroke-width="0.4"/>'
for i, b in enumerate(ba):
    x = 6 + (i + 0.5) / nb * (W5 - 10); yy = mid - b["beta"] / mx * (H5 - 22) / 2
    s5b += f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{mid}" y2="{yy:.1f}" stroke="{TEAL}" stroke-width="0.55"/>'
s5b += "</svg>"

# ---------- Step 6: ranking (top 30) and shrinkage
rk = d["rank"][:30]
W6, H6 = 300, 128; lo, hi = -0.1, 1.5
sy6 = lambda v: 12 + (hi - v) / (hi - lo) * (H6 - 34)
s6 = f'<svg viewBox="0 0 {W6} {H6}" width="100%" font-family="Inter">'
nsel = 24
s6 += f'<rect x="28" y="8" width="{(nsel)/30*(W6-34):.1f}" height="{H6-30}" fill="#FFF7ED"/>'
for t in [0, 0.5, 1.0, 1.5]:
    s6 += f'<line x1="28" x2="{W6-4}" y1="{sy6(t):.1f}" y2="{sy6(t):.1f}" stroke="#E2E8F0" stroke-width="0.5"/><text x="24" y="{sy6(t)+2.5:.1f}" text-anchor="end" font-size="8.0" fill="{MUTE}">{f(t,1)}</text>'
for i, r in enumerate(rk):
    x = 28 + (i + 0.5) / 30 * (W6 - 34); c = ORANGE if r["tested"] else INK
    s6 += f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{sy6(r["gebv"]-1.96*r["se"]):.1f}" y2="{sy6(r["gebv"]+1.96*r["se"]):.1f}" stroke="{c}" stroke-width="1" opacity="0.45"/>'
    s6 += f'<circle cx="{x:.1f}" cy="{sy6(r["gebv"]):.1f}" r="2.4" fill="{c if r["tested"] else "#fff"}" stroke="{c}" stroke-width="1.1"/>'
s6 += f'<text x="{28+nsel/30*(W6-34)/2:.1f}" y="20" text-anchor="middle" font-size="9.0" font-weight="800" fill="{ORANGE}">selected as parents (top 10%)</text>'
s6 += f'<text x="10" y="{H6/2-10}" text-anchor="middle" font-size="8.7" fill="{SLATE}" transform="rotate(-90 10 {H6/2-10})">GEBV (t ha⁻¹) ± 1.96 SE</text>'
s6 += f'<text x="{W6/2+12}" y="{H6-6}" text-anchor="middle" font-size="8.7" fill="{SLATE}">Rank (top 30 of {d["n"]})</text>'
s6 += "</svg>"

sh = d["shrink"]; h2 = d["h2"]
W7, H7 = 210, 150; lo7, hi7 = -1.9, 1.9
sx7 = lambda v: 30 + (v - lo7) / (hi7 - lo7) * (W7 - 38)
sy7 = lambda v: 8 + (hi7 - v) / (hi7 - lo7) * (H7 - 30)
s7 = f'<svg viewBox="0 0 {W7} {H7}" width="100%" font-family="Inter">'
s7 += f'<rect x="30" y="8" width="{W7-38}" height="{H7-30}" fill="#F8FAFC" stroke="#E2E8F0"/>'
for t in [-1, 0, 1]:
    s7 += f'<text x="26" y="{sy7(t)+2.5:.1f}" text-anchor="end" font-size="8.0" fill="{MUTE}">{f(t,0)}</text><text x="{sx7(t):.1f}" y="{H7-13}" text-anchor="middle" font-size="8.0" fill="{MUTE}">{f(t,0)}</text>'
s7 += f'<line x1="{sx7(lo7)}" y1="{sy7(lo7)}" x2="{sx7(hi7)}" y2="{sy7(hi7)}" stroke="{GRAY}" stroke-dasharray="3 2" stroke-width="0.9"/>'
s7 += f'<text x="{sx7(-1.75):.1f}" y="{sy7(-1.35):.1f}" font-size="8.3" fill="{MUTE}">1:1 = no shrinkage</text>'
for r in sh:
    s7 += f'<circle cx="{sx7(r["dev"]):.1f}" cy="{sy7(r["gebv"]):.1f}" r="2" fill="{TEAL}" opacity="0.6"/>'
s7 += f'<line x1="{sx7(lo7)}" y1="{sy7(lo7*h2):.1f}" x2="{sx7(hi7)}" y2="{sy7(hi7*h2):.1f}" stroke="{ORANGE}" stroke-width="1.8"/>'
s7 += f'<text x="{sx7(0.05):.1f}" y="{sy7(-0.75):.1f}" font-size="9.2" font-weight="800" fill="{ORANGE}">slope = h² = {h2:.2f}</text>'
ex = d["ex"]; exd = ex["ybar"] - ex["mean_y"]
s7 += f'<circle cx="{sx7(exd):.1f}" cy="{sy7(ex["gebv"]):.1f}" r="3.6" fill="none" stroke="{INK}" stroke-width="1.2"/>'
s7 += f'<text x="{sx7(exd)-5:.1f}" y="{sy7(ex["gebv"])-6:.1f}" text-anchor="end" font-size="8.5" font-weight="800" fill="{INK}">{ex["line"]}</text>'
s7 += f'<text x="{(W7+30)/2}" y="{H7-2}" text-anchor="middle" font-size="8.7" fill="{SLATE}">Adjusted mean − trial mean (t ha⁻¹)</text>'
s7 += f'<text x="9" y="{H7/2-6}" text-anchor="middle" font-size="8.7" fill="{SLATE}" transform="rotate(-90 9 {H7/2-6})">GEBV (t ha⁻¹)</text>'
s7 += "</svg>"

lam = d["se2"] / d["sg2"]
rep = dict(S1=s1, S2=s2, S3P=s3p, S3G=s3g, S4F=s4f, S4T=s4t, S5B=s5b, S6=s6, S7=s7,
           Y=vec(Y), Y2=vec(Y), XM=mat(X, "xm sm", color=xcol), BETA=vec([b * 1000 for b in B], 1),
           BETAH=vec([b * 1000 for b in B], 1), GM=mat(G, "gm sm", fmt=lambda v: f(v), color=gcol),
           GRR=rowvec(GR), GGB=rowvec(GV), LAM=f"{lam:.2f}", SG2=f"{d['sg2']:.3f}", SE2=f"{d['se2']:.3f}", H2=f"{h2:.2f}",
           MU=f"{d['mu']:.2f}", N=str(d["n"]), NT=str(d["n_train"]), NSNP=str(d["n_snp"]), ACC=f"{d['acc_untested']:.2f}",
           CORBB=f"{d['cor_bb']:.3f}", EXL=ex["line"], EXY=f"{ex['ybar']:.2f}", EXM=f"{ex['mean_y']:.2f}", EXD=f"{exd:.2f}",
           EXH=f"{exd*h2:.2f}", EXG=f"{ex['gebv']:.2f}", LROWS="".join(f"<div>{l}</div>" for l in L))
for tpl, out in (("fig1a_template.html", "Fig1_composite.html"), ("fig1b_template.html", "Fig2_composite.html")):
    html = open(os.path.join(here, tpl)).read()
    for k, v in rep.items(): html = html.replace("{{" + k + "}}", v)
    open(os.path.join(here, out), "w").write(html)
print("ok")
