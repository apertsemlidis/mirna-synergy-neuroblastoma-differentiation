"""JBS figure style helpers. See jbs.mplstyle for typography/axes defaults."""

from pathlib import Path
import matplotlib.pyplot as plt

_MM = 25.4
WIDTHS_IN = {
    "single": 85 / _MM,
    "1.5": 114 / _MM,
    "double": 170 / _MM,
}
MAX_HEIGHT_IN = 225 / _MM

# Wong (2011) colorblind-safe palette; keys match existing usage in this project.
PALETTE = {
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "green": "#009E73",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "black": "#000000",
}

STYLE_PATH = Path(__file__).with_name("jbs.mplstyle")


def use_style():
    plt.style.use(str(STYLE_PATH))


def setup_figure(width="single", height=None, nrows=1, ncols=1, **kwargs):
    use_style()
    w = WIDTHS_IN[width] if isinstance(width, str) else float(width)
    if height is None:
        height = min(w / 1.618, MAX_HEIGHT_IN)
    else:
        height = min(float(height), MAX_HEIGHT_IN)
    return plt.subplots(nrows, ncols, figsize=(w, height), **kwargs)


def save(fig, path, formats=("pdf", "png"), **kwargs):
    path = Path(path)
    stem = path.with_suffix("")
    for fmt in formats:
        fig.savefig(stem.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight", **kwargs)
