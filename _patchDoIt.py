#!/usr/bin/env python3
"""
Patch doIt.m:
  1. Remove all vessel-circle overlay plot lines
  2. Insert font-size setters (tick=12, text=16) before every save block
"""
import re, os

PATH = '/scratch/bass/projects/multiVencISMRM2026/doIt.m'

with open(PATH) as f:
    lines = f.readlines()

out = []
for line in lines:
    s = line.rstrip('\n')

    # ── 1. Remove circle overlay lines ─────────────────────────────────

    # "hold(ax,'on'); plot(... cVes ...)" on same line
    if re.search(r"hold\(ax\S*,'on'\).*plot\(.*cVes", s):
        continue
    # "hold(ax,'on'); plot(... cx, cy ...)" on same line
    if re.search(r"hold\(ax\S*,'on'\).*plot\(.*cx, cy", s):
        continue
    # "hold(axMask,'on'); plot(axMask, ID/2*cos..." on same line
    if re.search(r"hold\(axMask\S*,'on'\).*plot\(axMask", s):
        continue
    # standalone plot lines with ID/2 * cos(theta) or OD/2 * cos(theta)
    if re.search(r"^\s*plot\(ax\S*,\s*[IO]D/2\s*\*\s*cos\(theta", s):
        continue
    # standalone plot lines with axMask and cos(theta3)
    if re.search(r"^\s*plot\(axMask\S*,\s*[IO]D/2\s*\*cos\(theta3\)", s):
        continue

    # Remove "; hold on; plot(cx,cy,'w--','LineWidth',1)" from end of line
    s = re.sub(r";\s*hold on;\s*plot\(cx,cy,'w--','LineWidth',1\)", '', s)

    # ── 2. Font-size injection before save blocks ───────────────────────
    # Detect the figure handle for the current save block
    out.append(s + '\n')

# Second pass: inject font-size lines before every save block
# Pattern: "if saveThis || ~exist(fullfile(secXfig,'...'),'file')"
# and also "if saveThis || ~exist(fullfile(subFigDir,..."
# We inject BEFORE those lines, using the most recently seen figure handle.

SAVE_PAT  = re.compile(r"^\s*if saveThis \|\|")
FIG_PATS  = [
    (re.compile(r"\b(fSim)\s*=\s*figure"),     'fSim'),
    (re.compile(r"\b(fComb)\s*=\s*figure"),    'fComb'),
    (re.compile(r"\b(fMask)\s*=\s*figure"),    'fMask'),
    (re.compile(r"\b(f_ivs)\s*=\s*figure"),    'f_ivs'),
    (re.compile(r"\b(fVelSpec)\s*=\s*figure"), 'fVelSpec'),
    (re.compile(r"\b(fVelSpecInflow)\s*=\s*figure"), 'fVelSpecInflow'),
    (re.compile(r"\b(hFroi)\s*=\s*figure"),    'hFroi'),
    (re.compile(r"\b(hFv)\s*=\s*figure"),      'hFv'),
    (re.compile(r"\b(hF)\s*=\s*figure"),       'hF'),
    # generic "f = figure" — must come last
    (re.compile(r"^\s*\bf\b\s*=\s*figure"),    'f'),
]

current_fig = None
result = []
already_injected = set()   # track injection points to avoid duplicates

for i, line in enumerate(out):
    # update current figure handle
    for pat, handle in FIG_PATS:
        if pat.search(line):
            current_fig = handle
            break

    # inject before save block (once per consecutive save block)
    if SAVE_PAT.match(line) and current_fig and i not in already_injected:
        # look back to check we haven't just injected
        prev = result[-1].strip() if result else ''
        if not prev.startswith('set(findall('):
            indent = re.match(r'(\s*)', line).group(1)
            result.append(f"{indent}set(findall({current_fig},'Type','axes'),'FontSize',12);\n")
            result.append(f"{indent}set(findall({current_fig},'Type','text'),'FontSize',16);\n")
        already_injected.add(i)

    result.append(line)

with open(PATH, 'w') as f:
    f.writelines(result)

print(f"Done. {len(result)} lines written.")
