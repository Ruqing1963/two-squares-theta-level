"""Figure: sup and rms of (T - M)/sqrt(N) over N in (X/2, X] against mm', from sup_error_scan_x21.txt."""
import re
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "sup_error_scan_x21.txt"
rows = []
cur = None
for line in open(src, encoding="utf-8"):
    m = re.match(r"\(m, m'\) = \((\d+), (\d+)\), mm' = (\d+)", line)
    if m:
        cur = int(m.group(3))
        continue
    m = re.match(r"\s+all\s+n=\s*\d+\s+rms=\s*([\d.]+)\s+q99.9=\s*([\d.]+)\s+sup=\s*([\d.]+)", line)
    if m and cur is not None:
        rows.append((cur, float(m.group(1)), float(m.group(2)), float(m.group(3))))
        cur = None
rows = [r for r in rows if r[3] > 0]          # drop the exact cases (levels 4, 12)
rows.sort()
mm = [r[0] for r in rows]

fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.loglog(mm, [r[3] for r in rows], "o", color="#2a78d6", label=r"$\sup_N |T-M|/\sqrt{N}$")
ax.loglog(mm, [r[2] for r in rows], "s", color="#eb6834", label="99.9% quantile")
ax.loglog(mm, [r[1] for r in rows], "^", color="#1baf7a", label="root mean square")
x0, x1 = min(mm), max(mm)
ax.loglog([x0, x1], [1.76 * (x0 / 5) ** -0.25, 1.76 * (x1 / 5) ** -0.25], "--", color="#1baf7a", lw=1,
          label=r"$\propto (mm')^{-1/4}$")
ax.loglog([x0, x1], [14.7 * (x0 / 5) ** 0.25, 14.7 * (x1 / 5) ** 0.25], ":", color="#2a78d6", lw=1,
          label=r"$\propto (mm')^{1/4}$ (Theorem 2.3)")
ax.set_xlabel(r"$mm'$")
ax.set_ylabel(r"$(T-M)/\sqrt{N}$,  $N \in (2^{20}, 2^{21}]$")
ax.grid(True, which="both", alpha=0.25)
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout()
fig.savefig("sup_error.pdf")
fig.savefig("sup_error.png", dpi=150)
print(f"{len(rows)} points; wrote sup_error.pdf/.png")
