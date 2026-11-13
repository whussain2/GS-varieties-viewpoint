# Builds Fig2_composite.html (Box 2). Panel C uses ../code/Fig3C_locus_frequencies.csv from the R simulation.
import csv, os, random
here = os.path.dirname(os.path.abspath(__file__))
TEAL, TEAL_L, ORANGE, SLATE, MUTE, INK, GRAY = "#0F766E", "#99F6E4", "#C2410C", "#334155", "#64748B", "#0F172A", "#94A3B8"

# loci positions along a 70-unit chromosome; chr1: tight pair (l1,l2) + loose l3; chr2: three spread loci
POS = {1: [8, 17, 56], 2: [9, 35, 58]}
def homolog(x, y, alleles, chrn, hl=None):
    s = f'<rect x="{x}" y="{y}" width="66" height="7" rx="3.5" fill="#E2E8F0" stroke="#CBD5E1" stroke-width="0.6"/>'
    if hl:  # highlight a segment (recombinant part)
        a, b = hl; s += f'<rect x="{x+a}" y="{y}" width="{b-a}" height="7" rx="3.5" fill="#FED7AA" opacity="0.8"/>'
    for k, a in enumerate(alleles):
        cx = x + POS[chrn][k]; cy = y + 3.5
        if a == "+":
            s += f'<circle cx="{cx}" cy="{cy}" r="3.8" fill="{TEAL}"/><text x="{cx}" y="{cy+2.6}" text-anchor="middle" font-size="7.8" font-weight="800" fill="#fff">+</text>'
        else:
            s += f'<circle cx="{cx}" cy="{cy}" r="3.8" fill="#fff" stroke="{GRAY}" stroke-width="1"/><text x="{cx}" y="{cy+2.4}" text-anchor="middle" font-size="8.4" font-weight="800" fill="{GRAY}">−</text>'
    return s

def plant(x, y, name, g, sel=False, w=158):
    c1a, c1b, c2a, c2b = g
    n = sum(s.count("+") for s in g)
    s = f'<g transform="translate({x},{y})">'
    s += f'<rect x="0" y="0" width="{w}" height="40" rx="6" fill="{"#FFF7ED" if sel else "#fff"}" stroke="{ORANGE if sel else "#E2E8F0"}" stroke-width="{1.4 if sel else 0.8}"/>'
    s += f'<text x="6" y="10" font-size="8.3" font-weight="800" fill="{INK}">{name}</text>'
    s += f'<text x="{w-6}" y="10" text-anchor="end" font-size="8.3" font-weight="800" fill="{ORANGE if sel else SLATE}">{n}/12 +</text>'
    s += homolog(6, 15, c1a, 1) + homolog(6, 26, c1b, 1) + homolog(84, 15, c2a, 2) + homolog(84, 26, c2b, 2)
    return s + "</g>"

pop = [("Plant 1", ["+--", "-+-", "+-+", "-++"]), ("Plant 2", ["+-+", "-+-", "++-", "-++"]),
       ("Plant 3", ["---", "-+-", "--+", "+--"]), ("Plant 4", ["++-", "--+", "-++", "++-"]),
       ("Plant 5", ["-+-", "--+", "---", "++-"]), ("Plant 6", ["+--", "--+", "-+-", "+-+"])]
# ---------------- A
A = '<svg viewBox="0 0 500 128" width="100%" font-family="Inter">'
for col in range(3):
    cx = 2 + col * 166
    A += f'<text x="{cx+39}" y="8" text-anchor="middle" font-size="7.4" fill="{MUTE}">chromosome 1</text><text x="{cx+117}" y="8" text-anchor="middle" font-size="7.4" fill="{MUTE}">chromosome 2</text>'
for i, (nm, g) in enumerate(pop):
    A += plant(2 + (i % 3) * 166, 12 + (i // 3) * 46, nm, g, sel=nm in ("Plant 2", "Plant 4"))
A += f'<g transform="translate(2,110)">'
A += f'<circle cx="5" cy="6" r="3.8" fill="{TEAL}"/><text x="5" y="8.6" text-anchor="middle" font-size="7.8" font-weight="800" fill="#fff">+</text><text x="13" y="8.5" font-size="8.3" fill="{SLATE}">allele that increases yield</text>'
A += f'<circle cx="118" cy="6" r="3.8" fill="#fff" stroke="{GRAY}"/><text x="118" y="8.4" text-anchor="middle" font-size="8.4" font-weight="800" fill="{GRAY}">−</text><text x="126" y="8.5" font-size="8.3" fill="{SLATE}">allele that decreases yield</text>'
A += f'<rect x="232" y="1" width="14" height="10" rx="3" fill="#FFF7ED" stroke="{ORANGE}" stroke-width="1.2"/><text x="250" y="8.5" font-size="8.3" fill="{ORANGE}" font-weight="700">selected as parents</text>'
A += f'<text x="340" y="8.5" font-size="8.3" fill="{SLATE}"><tspan font-weight="800" fill="{INK}">Loci 1–2:</tspan> tightly linked</text>'
A += '</g></svg>'

# ---------------- B : cross, meiosis, recombinant progeny (left to right)
def pair(x, y, a, b, chrn, hla=None, hlb=None): return homolog(x, y, a, chrn, hla) + homolog(x, y + 10, b, chrn, hlb)
B = '<svg viewBox="0 0 690 92" width="100%" font-family="Inter">'
B += f'<defs><marker id="ar" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="{GRAY}"/></marker></defs>'
B += f'<rect x="2" y="2" width="172" height="86" rx="7" fill="#FFF7ED" stroke="{ORANGE}" stroke-width="1"/>'
B += f'<text x="10" y="13" font-size="8.5" font-weight="800" fill="{INK}">Crossing block</text>'
B += f'<text x="10" y="26" font-size="8.1" font-weight="800" fill="{ORANGE}">Plant 2</text><text x="166" y="26" text-anchor="end" font-size="8.1" font-weight="700" fill="{ORANGE}">7/12 +</text>' + pair(10, 30, "+-+", "-+-", 1) + pair(92, 30, "++-", "-++", 2)
B += f'<text x="88" y="56" text-anchor="middle" font-size="12.3" font-weight="700" fill="{INK}">×</text>'
B += f'<text x="10" y="62" font-size="8.1" font-weight="800" fill="{ORANGE}">Plant 4</text><text x="166" y="62" text-anchor="end" font-size="8.1" font-weight="700" fill="{ORANGE}">7/12 +</text>' + pair(10, 66, "++-", "--+", 1) + pair(92, 66, "-++", "++-", 2)
B += f'<path d="M180 45 H214" stroke="{GRAY}" stroke-width="1.4" marker-end="url(#ar)"/><text x="197" y="38" text-anchor="middle" font-size="7.6" fill="{MUTE}">meiosis</text>'
B += f'<rect x="220" y="2" width="172" height="86" rx="7" fill="#fff" stroke="#E2E8F0"/>'
B += f'<text x="228" y="13" font-size="8.5" font-weight="800" fill="{INK}">Gametes</text><text x="384" y="13" text-anchor="end" font-size="7.4" fill="{MUTE}">shaded: crossover segment</text>'
B += f'<text x="228" y="30" font-size="7.8" fill="{MUTE}">from Plant 2</text>' + homolog(228, 34, "+-+", 1) + homolog(310, 34, "+++", 2, (23, 66))
B += f'<text x="228" y="62" font-size="7.8" fill="{MUTE}">from Plant 4</text>' + homolog(228, 66, "+++", 1, (36, 66)) + homolog(310, 66, "+++", 2, (0, 23))
B += f'<path d="M398 45 H432" stroke="{GRAY}" stroke-width="1.4" marker-end="url(#ar)"/><text x="415" y="38" text-anchor="middle" font-size="7.6" fill="{MUTE}">fertilization</text>'
B += f'<rect x="438" y="14" width="168" height="54" rx="7" fill="#F0FDFA" stroke="{TEAL}" stroke-width="1.4"/>'
B += f'<text x="446" y="27" font-size="8.5" font-weight="800" fill="{INK}">Best progeny</text><text x="598" y="27" text-anchor="end" font-size="8.5" font-weight="800" fill="{TEAL}">11/12 +</text>'
B += pair(446, 34, "+-+", "+++", 1) + pair(528, 34, "+++", "+++", 2)
B += f'<text x="612" y="30" font-size="8.1" fill="{SLATE}">better than</text><text x="612" y="40" font-size="8.1" fill="{SLATE}">either parent</text>'
B += f'<text x="612" y="56" font-size="8.5" font-weight="800" fill="{ORANGE}">≈1 in 1,200</text><text x="612" y="66" font-size="8.1" fill="{SLATE}">progeny</text>'
B += f'<text x="438" y="82" font-size="7.4" fill="{SLATE}">The − tightly linked to a + (chromosome 1) stays.</text>'
B += "</svg>"

# ---------------- C : allele frequencies over cycles (simulation)
rows = list(csv.DictReader(open(os.path.join(here, "../code/Fig3C_locus_frequencies.csv"))))
loci = {}
for r in rows: loci.setdefault(int(r["locus"]), []).append((int(r["cycle"]), float(r["freq"]), float(r["abs_effect"])))
nc = 12; W, H = 250, 250; x0, y0, x1, y1 = 30, 8, 244, 222
sx = lambda c: x0 + c / nc * (x1 - x0); sy = lambda f: y1 - f * (y1 - y0)
C = f'<svg viewBox="0 0 {W} {H+20}" width="100%" font-family="Inter">'
C += f'<rect x="{x0}" y="{y0}" width="{x1-x0}" height="{y1-y0}" fill="#F8FAFC" stroke="#E2E8F0"/>'
for t in [0, 0.25, 0.5, 0.75, 1]:
    C += f'<line x1="{x0}" x2="{x1}" y1="{sy(t)}" y2="{sy(t)}" stroke="#E2E8F0" stroke-width="0.5"/><text x="{x0-3}" y="{sy(t)+2.5}" text-anchor="end" font-size="7.6" fill="{MUTE}">{int(t*100)}%</text>'
for c in range(0, 13, 4): C += f'<text x="{sx(c)}" y="{y1+10}" text-anchor="middle" font-size="7.6" fill="{MUTE}">{c}</text>'
random.seed(4); pick = random.sample(sorted(loci), 36)
emax = max(v[0][2] for v in loci.values())
for k in sorted(pick, key=lambda k: loci[k][0][2]):
    tr = sorted(loci[k]); e = tr[0][2] / emax
    pts = " ".join(f"{sx(c):.1f},{sy(f):.1f}" for c, f, _ in tr)
    C += f'<polyline points="{pts}" fill="none" stroke="{TEAL}" stroke-width="{0.5+1.3*e:.2f}" opacity="{0.18+0.6*e:.2f}"/>'
mean = [sum(dict((c, f) for c, f, _ in loci[k])[c] for k in loci) / len(loci) for c in range(nc + 1)]
C += f'<polyline points="{" ".join(f"{sx(c):.1f},{sy(m):.1f}" for c, m in enumerate(mean))}" fill="none" stroke="{ORANGE}" stroke-width="2.4"/>'
C += f'<text x="{sx(12)-2}" y="{sy(mean[-1])+12}" text-anchor="end" font-size="8.7" font-weight="800" fill="{ORANGE}" stroke="#fff" stroke-width="2.4" paint-order="stroke">mean of 200 loci: {mean[0]*100:.0f}% → {mean[-1]*100:.0f}%</text>'
fixed = sum(1 for k in loci if sorted(loci[k])[-1][1] == 1); lost = sum(1 for k in loci if sorted(loci[k])[-1][1] == 0)
C += f'<text x="{(x0+x1)/2}" y="{y1+21}" text-anchor="middle" font-size="8.1" fill="{SLATE}">Cycles of selection</text>'
C += f'<text x="8" y="{(y0+y1)/2}" text-anchor="middle" font-size="8.1" fill="{SLATE}" transform="rotate(-90 8 {(y0+y1)/2})">Frequency of the + allele</text>'
C += f'<text x="{x0}" y="{H+8}" font-size="7.6" fill="{MUTE}">Thin lines: single loci (thicker = larger effect).</text><text x="{x0}" y="{H+17}" font-size="7.6" fill="{MUTE}">After 12 cycles: {fixed} + alleles fixed, {lost} lost (mostly by drift).</text>'
C += "</svg>"

html = open(os.path.join(here, "fig2_template.html")).read()
for k, v in dict(A=A, B=B, C=C).items(): html = html.replace("{{" + k + "}}", v)
open(os.path.join(here, "Fig3_composite.html"), "w").write(html)
print("ok", fixed, lost, [round(m, 3) for m in mean])
