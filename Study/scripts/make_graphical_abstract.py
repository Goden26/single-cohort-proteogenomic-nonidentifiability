#!/usr/bin/env python
"""Graphical abstract for the JBI submission. Renders a single horizontal 5-panel schematic at the
JBI-required proportions (>= 531 x 1328 px h x w, readable at 5 x 13 cm). Outputs PNG + PDF.
Color logic: blue = batch-free / trustworthy, red = aliased false positives, green = surviving positive."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SUB=os.path.join(os.path.dirname(ROOT),"Submission package")
BLUE="#2c7fb8"; RED="#d7301f"; GREEN="#2ca25f"; GREY="#555555"; INK="#1a1a1a"

fig=plt.figure(figsize=(13.3/2.54, 5.3/2.54))           # ~13.3 x 5.3 cm
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,40); ax.axis("off")

def box(x,y,w,h,fc,ec,lw=1.2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.3,rounding_size=1.2",
                                fc=fc,ec=ec,lw=lw,mutation_aspect=0.5))
def arrow(x0,x1,y,label):
    ax.add_patch(FancyArrowPatch((x0,y),(x1,y),arrowstyle="-|>",mutation_scale=13,lw=2.0,color=INK))
    ax.text((x0+x1)/2,y+2.2,label,ha="center",va="bottom",fontsize=4.4,color=INK,style="italic")
def txt(x,y,s,**k): ax.text(x,y,s,**k)

txt(50,37.3,"Single-cohort-per-condition proteogenomic merges are non-identifiable for the cross-condition contrast",
    ha="center",va="center",fontsize=5.4,weight="bold",color=INK)

# PANEL 1 - two single cohorts
box(1.5,14,15,16,"#eaf3f8",BLUE)
txt(9,27.6,"single cohort\nper condition",ha="center",va="center",fontsize=4.6,weight="bold",color=BLUE)
box(4,21,11,3.0,"#ffffff",BLUE,0.8); txt(9.5,22.5,"HCC proteome",ha="center",va="center",fontsize=4.0,color=INK)
box(4,16.4,11,3.0,"#ffffff",BLUE,0.8); txt(9.5,17.9,"iCCA proteome",ha="center",va="center",fontsize=4.0,color=INK)
arrow(17,22,21.5,"merge")

# PANEL 2 - aliasing
box(22,14,16,16,"#fdecea",RED)
txt(30,27.6,"condition = batch",ha="center",va="center",fontsize=5.0,weight="bold",color=RED)
ax.add_patch(plt.Rectangle((24.5,20.2),11,2.4,fc=RED,ec="none",alpha=0.85))
ax.add_patch(plt.Rectangle((24.5,20.2),11,2.4,fc="none",ec=BLUE,lw=1.4,ls=(0,(3,2))))
txt(30,17.6,"perfectly aliased\n(non-identifiable)",ha="center",va="center",fontsize=4.2,color=INK)
arrow(38.5,44.5,22,"harmonize\n/ ComBat")

# PANEL 3 - false hits
box(45,14,16,16,"#fdecea",RED)
txt(53,27.6,"false hits",ha="center",va="center",fontsize=5.0,weight="bold",color=RED)
rng=np.random.default_rng(3)
ax.scatter(53+rng.normal(0,3.2,130),22+rng.normal(0,2.4,130),s=1.8,c=RED,alpha=0.55,linewidths=0)
txt(53,16.6,"63% liver / 76% kidney\nflagged significant",ha="center",va="center",fontsize=4.0,color=INK)
arrow(61.5,67,22,"validate")

# PANEL 4 - batch-free RNA anchor exposes them
box(67.5,14,15.5,16,"#eaf3f8",BLUE)
txt(75.2,28.0,"batch-free RNA anchor",ha="center",va="center",fontsize=4.4,weight="bold",color=BLUE)
axc=fig.add_axes([0.715,0.515,0.072,0.15]); axc.set_xticks([]); axc.set_yticks([])
for sp in axc.spines.values(): sp.set_color(GREY); sp.set_linewidth(0.6)
axc.axhline(0,color=GREY,lw=0.4); axc.axvline(0,color=GREY,lw=0.4)
axc.scatter(rng.normal(0,1,110),rng.normal(0,1,110),s=1.4,c=GREY,alpha=0.5,linewidths=0)
axc.set_xlim(-3,3); axc.set_ylim(-3,3)
txt(75.2,17.0,"none validate (rho = -0.15)",ha="center",va="center",fontsize=3.9,color=RED)
txt(75.2,15.2,"except LAT1 (trial-anchored)",ha="center",va="center",fontsize=3.7,color=GREEN,weight="bold")

# PANEL 5 - calibration curve
box(83.5,14,15,16,"#eef7f1",GREEN)
txt(91,28.0,"anchor check calibrated",ha="center",va="center",fontsize=4.4,weight="bold",color=GREEN)
axk=fig.add_axes([0.862,0.515,0.105,0.15])
a=np.linspace(0,1,11)
axk.plot(a,0.97-0.12*a,"-",color=BLUE,lw=1.3)
axk.plot(a,0.06+0.87*a,"-",color=RED,lw=1.3)
axk.plot(a,0.42-0.18*a,"--",color=GREY,lw=1.0)
axk.set_xlim(0,1); axk.set_ylim(0,1); axk.set_xticks([]); axk.set_yticks([])
for sp in axk.spines.values(): sp.set_color(GREY); sp.set_linewidth(0.6)
txt(91,17.0,"recall stays high (blue);",ha="center",va="center",fontsize=3.7,color=INK)
txt(91,15.4,"false discovery to 0.93 (red)",ha="center",va="center",fontsize=3.7,color=INK)

# footer ribbon
ax.add_patch(plt.Rectangle((1.5,8.5),97,3.4,fc="#f0f0f0",ec=GREY,lw=0.6))
txt(50,10.2,"Reproducible benchmark: liver (powered) + kidney (corroborating), all open data and code",
    ha="center",va="center",fontsize=4.4,color=INK,weight="bold")

for ext in ["png","pdf"]:
    fig.savefig(os.path.join(SUB,f"graphical_abstract.{ext}"),dpi=300)   # exact figsize, preserves ~2.5:1
print("saved graphical_abstract.png and .pdf to", SUB)
