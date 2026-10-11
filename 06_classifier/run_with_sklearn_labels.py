"""Launcher that replaces the internal 'classif.xxx' labels with scikit-learn class names before figures are saved.
Usage: python run_with_sklearn_labels.py <script path> [script arguments...]"""
import re, runpy, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "06_classifier")
from fig6_model_labels import MODEL_LABELS
PAT = re.compile(r"classif\.[a-z_]+")
def disp(s):
    return PAT.sub(lambda m: MODEL_LABELS.get(m.group(0), {"display": m.group(0)})["display"], s)
_orig = mf.Figure.savefig
def _savefig(self, *a, **k):
    self.canvas.draw()
    for ax in self.axes:
        for getter, setter in ((ax.get_xticklabels, ax.set_xticklabels), (ax.get_yticklabels, ax.set_yticklabels)):
            labs = [t.get_text() for t in getter()]
            if any("classif." in l for l in labs):
                setter([disp(l) for l in labs])
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                t.set_text(disp(t.get_text()))
        for t in ax.texts:
            t.set_text(disp(t.get_text()))
    for t in self.texts:
        t.set_text(disp(t.get_text()))
    return _orig(self, *a, **k)
mf.Figure.savefig = _savefig
script = sys.argv[1]
sys.argv = [script] + sys.argv[2:]
runpy.run_path(script, run_name="__main__")
