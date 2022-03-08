import os
from flask import (
    Flask,
    make_response,
    request,
    redirect,
    url_for)
from jinja2 import Environment, FileSystemLoader, Template
import utils
import legendrian_links as ll

LINKS = utils.LINKS
DEFAULT_LINK_KEY = 'm10_132_nv'

# Static resources
STATIC_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'static')
JINJA_ENV = Environment(loader=FileSystemLoader(STATIC_PATH))
INDEX_TEMPLATE = JINJA_ENV.get_template('index.html')

# App resources
APP = Flask(__name__)

KNOT_COLORS = [
    (255, 0, 0),
    (128, 0, 128),
    (0, 0, 255),
    (0, 128, 128),
    (0, 255, 0),
    (128, 128, 0)
]
KNOT_ORIENTATIONS_TO_ARROW = {
    'l': '<',
    'r': '>'
}


def get_template_context(pd):
    chords = [
        {
            "string": str(chord),
            "from_knot": chord.bottom_line_segment.knot_label,
            "to_knot": chord.top_line_segment.knot_label
        }
        for chord in pd.chords]
    dgas = []
    if not pd.lazy_lch:
        dgas.append(get_dga_context(pd.lch_dga, name="LCH"))
    if not pd.lazy_rsft:
        dgas.append(get_dga_context(pd.rsft_dga, name="RSFT"))
    template_context = {
        'front_crossings': ','.join([str(c) for c in pd.front_crossings]),
        'svg_context': get_diagram_context(pd),
        'chords': chords,
        'dgas': dgas
    }
    return template_context


def get_diagram_context(pd, increment=50, pad=10):
    knots = pd.knots
    for k in range(len(knots)):
        knots[k]["label"] = k
        knots[k]["rgb"] = KNOT_COLORS[k % len(KNOT_COLORS)]
    n_strands = pd.n_strands
    n_segments = pd.max_x_right
    height = increment * n_strands + pad
    width = (increment * n_segments) + pad + increment
    line_segments = pd.get_line_segment_array()
    output = {
        "knots": knots,
        "pad": pad,
        "increment": increment,
        "height": height,
        "width": width,
        "lines": [
            {
                'start_xy': [
                    increment + pad + int(increment * ls['array'][0][0]), pad + int(increment * ls['array'][0][1])],
                'end_xy': [
                    increment + pad + int(increment * ls['array'][1][0]), pad + int(increment * ls['array'][1][1])],
                'rgb': KNOT_COLORS[ls['knot_label'] % len(KNOT_COLORS)],
                'label': {
                    'x': increment + pad + int(increment * (ls['array'][0][0] + ls['array'][1][0]) / 2),
                    'y': pad + int(increment * (ls['array'][0][1] + ls['array'][1][1]) / 2),
                    'marker': KNOT_ORIENTATIONS_TO_ARROW[ls['orientation']] if ls['t_label'] else ''
                }
            }
            for ls in line_segments
        ],
        'x_labels': [{"label": x, "x": int(increment + increment * x + increment / 2)} for x in range(n_segments)],
        'y_labels': [{"label": y, "y": int(pad + increment * y + increment / 2)} for y in range(n_strands - 1)],
    }
    return output


def get_dga_context(dga, name):
    output = dict()
    output["name"] = name
    output["grading_mod"] = dga.grading_mod
    output["coeff_mod"] = dga.coeff_mod
    generators = [
        {
            "symbol": g,
            "name": str(g),
            "grading": dga.gradings[g],
            "del": str(dga.differentials[g])
        }
        for g in dga.symbols
    ]
    generators = sorted(generators, key=lambda g: g["name"])
    output["generators"] = generators
    output["n_augs"] = len(dga.augmentations)
    output["has_augs"] = len(dga.augmentations) > 0
    deg_0_gens = [g for g in generators if g["grading"] == 0]
    output["deg_0_gens"] = deg_0_gens
    output["augs"] = [{g["name"]: aug[g["symbol"]] for g in deg_0_gens} for aug in dga.augmentations]
    output["bilin_polys"] = dga.bilin_polys
    return output


@APP.route('/')
def home():
    n_strands = request.args.get('n_strands')
    if n_strands is not None:
        n_strands = int(n_strands)
        crossings = request.args.get('crossings')
        if crossings is None:
            crossings = []
        else:
            crossings = [int(x) for x in crossings.split(",")]
        lazy_lch = True
        lazy_rsft = True
        auto_dga = request.args.get('auto_dgas')
        if auto_dga is not None:
            auto_dgas = auto_dga.lower().split(',')
            if 'lch' in auto_dgas:
                lazy_lch = False
            if 'rsft' in auto_dgas:
                lazy_rsft = False
        n_copy = request.args.get('n_copy')
        if n_copy is not None:
            n_copy = int(n_copy)
        else:
            n_copy = 1
        lazy_disks = False
        lazy_disks_flag = request.args.get('lazy_disks')
        if lazy_disks_flag is not None:
            lazy_disks = lazy_disks_flag.lower() == 'true'
        pd = ll.PlatDiagram(
            n_strands=n_strands,
            front_crossings=crossings,
            n_copy=n_copy,
            lazy_disks=lazy_disks,
            lazy_lch=lazy_lch,
            lazy_rsft=lazy_rsft,
        )
    else:
        data = LINKS[DEFAULT_LINK_KEY].copy()
        data.pop('comment')
        data['lazy_lch'] = True
        data['lazy_rsft'] = False
        pd = ll.PlatDiagram(**data)
    template_context = get_template_context(pd=pd)
    context = INDEX_TEMPLATE.render(**template_context)
    return make_response(context)


if __name__ == "__main__":
    APP.run()
