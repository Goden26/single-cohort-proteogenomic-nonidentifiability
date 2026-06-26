#!/usr/bin/env python
"""Rebuild the JBI submission deliverables from the markdown sources (portable; run from anywhere).
Produces the separate files Elsevier/JBI wants plus a combined convenience PDF:
  - manuscript.docx        (figures embedded under their legends; Times New Roman reference doc)
  - cover_letter.docx, title_page.docx, declarations.docx, suggested_reviewers.docx
  - graphical_abstract.png/.pdf are produced separately by make_graphical_abstract.py
  - submission_package.pdf  (title + cover + graphical abstract + manuscript + declarations)
Requires pandoc and tectonic on PATH."""
import os, re, subprocess, tempfile, glob
HERE=os.path.dirname(os.path.abspath(__file__)); STUDY=os.path.dirname(HERE); BASE=os.path.dirname(STUDY)
SUB=os.path.join(BASE,"Submission package"); FIG=os.path.join(STUDY,"figures"); REF=os.path.join(BASE,"Materials","_tnr_ref.docx")
GA=os.path.join(SUB,"graphical_abstract.png")
figmap={n:sorted(glob.glob(os.path.join(FIG,f"fig{n}_*.png")))[0] for n in range(1,7)}

def with_figures(md):
    def repl(m):
        n=int(m.group(1)); img=figmap.get(n)
        return m.group(0)+(f"\n\n![]({img}){{width=85%}}\n" if img else "")
    return re.sub(r"\*\*Figure (\d)\.[^\n]*\*\*[^\n]*", repl, md)

man=open(os.path.join(SUB,"manuscript.md")).read()
with tempfile.NamedTemporaryFile("w",suffix=".md",delete=False,dir=SUB) as tf:
    tf.write(with_figures(man)); man_fig=tf.name

# 1) manuscript.docx (with figures)
subprocess.run(["pandoc",man_fig,"-o",os.path.join(SUB,"manuscript.docx"),
                "--reference-doc",REF,"--resource-path",FIG],check=True)
print("built manuscript.docx")

# 2) individual submission docx files
for stem in ["cover_letter","title_page","declarations","suggested_reviewers"]:
    src=os.path.join(SUB,f"{stem}.md")
    if os.path.exists(src):
        subprocess.run(["pandoc",src,"-o",os.path.join(SUB,f"{stem}.docx"),"--reference-doc",REF],check=True)
        print(f"built {stem}.docx")

# 3) combined submission_package.pdf (title + cover + graphical abstract + manuscript + declarations)
with tempfile.NamedTemporaryFile("w",suffix=".md",delete=False,dir=SUB) as gf:
    ga_md = f"# Graphical abstract\n\n![]({GA})\n\n" if os.path.exists(GA) else ""
    gf.write(ga_md); ga_file=gf.name
parts=[os.path.join(SUB,"title_page.md"),os.path.join(SUB,"cover_letter.md"),ga_file,man_fig,os.path.join(SUB,"declarations.md")]
parts=[p for p in parts if os.path.exists(p)]
try:
    subprocess.run(["pandoc",*parts,"-o",os.path.join(SUB,"submission_package.pdf"),
                    "--pdf-engine=tectonic","--resource-path",f"{FIG}:{SUB}",
                    "-V","geometry:margin=1in","-V","fontsize=11pt"],check=True)
    print("built submission_package.pdf")
except subprocess.CalledProcessError as e:
    print("PDF build failed (docx files are still current):",e)
finally:
    os.unlink(man_fig); os.unlink(ga_file)
print("BUILD_DONE")
