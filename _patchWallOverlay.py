#!/usr/bin/env python3
"""
Patch doIt.m: add opaque gray wall-mask overlay to all velocity/phase imagesc maps,
fix velocity CLims to [-venc,venc], and insert grayRGB_mask definitions.
"""
import re

PATH = '/scratch/bass/projects/multiVencISMRM2026/doIt.m'
GRAY_DEF  = "grayRGB_mask = 0.45 * ones([size(maskWallLowMag), 3]);"
OVL_AX    = "hold({ax},'on'); h_ov=image({ax},PEpos,FEpos,grayRGB_mask); h_ov.AlphaData=double(maskWallLowMag);"

with open(PATH) as f:
    lines = f.readlines()

out = []
i = 0
while i < len(lines):
    line = lines[i]
    s = line.rstrip('\n')

    # ── Insert grayRGB_mask before the magPhaseMaps figure (sec5) ────────────
    if s.strip() == "% figures" and i+1 < len(lines) and "tiledlayout(f,4,6" in lines[i+1]:
        indent = re.match(r'(\s*)', s).group(1)
        out.append(f"{indent}{GRAY_DEF}\n")
        out.append(line)
        i += 1; continue

    # ── Insert grayRGB_mask before sim-setup in sec4 ─────────────────────────
    if s.strip() == "% --- Simulation: setup matched to joint fit ---":
        indent = re.match(r'(\s*)', s).group(1)
        out.append(f"{indent}{GRAY_DEF}\n")
        out.append(line)
        i += 1; continue

    # ── Insert grayRGB_mask before Mag-flow-on figure in sec9 ────────────────
    if s.strip() == "% Mag flow on":
        indent = re.match(r'(\s*)', s).group(1)
        out.append(f"{indent}{GRAY_DEF}\n")
        out.append(line)
        i += 1; continue

    # ── sec5 magPhaseMaps: phase tiles — add overlay after title line ─────────
    if re.search(r"title\(ax\{end\},\s*['\[].*phase['\]]", s):
        indent = re.match(r'(\s*)', s).group(1)
        out.append(f"{indent}{OVL_AX.format(ax='ax{end}')}\n")
        out.append(line)
        i += 1; continue

    # ── sec4 fComb: phantom vel CLim ─────────────────────────────────────────
    if "imagesc(axComb{end}, PEpos, FEpos, im, [-1 1].* max(abs(im(:))))" in s:
        s = s.replace("[-1 1].* max(abs(im(:)))", "[-bestVenc bestVenc]")
        line = s + '\n'

    # ── sec4 fComb: add overlay after phantom-vel colormap set ───────────────
    if re.search(r"axComb\{end\}\.Colormap = blueBlackRed.*XColor.*'m'", s):
        out.append(line)
        indent = re.match(r'(\s*)', s).group(1)
        out.append(f"{indent}{OVL_AX.format(ax='axComb{end}')}\n")
        i += 1; continue

    # ── sec4 fMask: phantom vel CLim ─────────────────────────────────────────
    if "imagesc(axMask{end}, PEpos, FEpos, im, [-1 1].*max(abs(im(:))))" in s:
        s = s.replace("[-1 1].*max(abs(im(:)))", "[-bestVenc bestVenc]")
        line = s + '\n'

    # ── sec4 fMask: add overlay after phantom-vel colormap+title ─────────────
    if re.search(r"title\(axMask\{end\},\s*\[.phantom vel", s):
        indent = re.match(r'(\s*)', s).group(1)
        out.append(f"{indent}{OVL_AX.format(ax='axMask{end}')}\n")
        out.append(line)
        i += 1; continue

    # ── sec9 phase loops: overlay after Colormap=blueBlackRed (inside loop) ──
    # Detect: line is "    ax{end}.Colormap = blueBlackRed;" inside a phase loop
    # Heuristic: preceded by loop context (vencList) — use flag
    if re.search(r"ax\{end\}\.Colormap = blueBlackRed;$", s.rstrip()):
        # Check if this is inside a velocity/phase loop (look back for vencList loop)
        ctx = ''.join(l for l in out[-30:])
        if re.search(r'(phase2vel|angle\(mean\(data|angle\(mean\(dataNoFlow|CDvel|PDvel|Pvel)', ctx):
            out.append(line)
            indent = re.match(r'(\s*)', s).group(1)
            # For CDvel/PDvel/Pvel: also set per-venc CLim before overlay
            if re.search(r'(CDvel|PDvel|Pvel)', ctx):
                out.append(f"{indent}if isfinite(vencList(vencIdx)); set(ax{{end}},'CLim',[-vencList(vencIdx) vencList(vencIdx)]); else; set(ax{{end}},'CLim',[-9 9]); end\n")
            out.append(f"{indent}{OVL_AX.format(ax='ax{end}')}\n")
            i += 1; continue

    # ── sec9 CDvel/PDvel/Pvel: remove global set([ax{:}],'CLim',...) ─────────
    if re.search(r"set\(\[ax\{:\}\],'CLim',\[-[19] [19]\]", s):
        # Only remove if it's the CDvel/PDvel/Pvel post-loop setter
        # (leave others like mag-flow CLim setter intact)
        if 'cLim' not in s and 'max' not in s:
            i += 1; continue  # drop this line

    out.append(line)
    i += 1

with open(PATH, 'w') as f:
    f.writelines(out)

# Count injections
n_ovl  = sum(1 for l in out if 'h_ov=image' in l)
n_gray = sum(1 for l in out if 'grayRGB_mask = 0.45' in l)
print(f"grayRGB_mask defs: {n_gray},  overlay injections: {n_ovl}")
