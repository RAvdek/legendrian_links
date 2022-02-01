import os
from flask import (
    Flask,
    make_response,
    request,
    redirect,
    url_for)
from jinja2 import Template
import utils
import legendrian_links as ll

LINKS = utils.LINKS
DEFAULT_LINK_KEY = 'm10_132_nv'

# Static resources
STATIC_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'static')
HTML_PATH = os.path.join(STATIC_PATH, 'index.html')
with open(HTML_PATH) as f:
    INDEX_TEMPLATE = Template(f.read())

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


def get_template_context(pd, increment=50, pad=10):
    n_strands = pd.n_strands
    n_segments = len(pd.plat_segments)
    height = increment * pd.n_strands + pad
    width = (increment * n_segments) + pad + increment
    line_segments = pd.get_line_segment_array()
    disk_corners = pd.disk_corners
    knots = pd.knots
    for k in range(len(knots)):
        knots[k]["label"] = k
        knots[k]["rgb"] = KNOT_COLORS[k % len(KNOT_COLORS)]
    chords = [
        {
            "string": str(chord),
            "grading": str(pd.get_lch_generator_from_chord(chord).grading),
            "from_knot": chord.bottom_line_segment.knot_label,
            "to_knot": chord.top_line_segment.knot_label,
            "lch_del": str(pd.lch_del[pd.get_lch_generator_from_chord(chord).symbol].to_polynomial())
        }
        for chord in pd.chords]
    rsft_generators = [
        {
            "string": str(word),
            "grading": word.grading
        }
        for word in pd.rsft_generators
    ]
    template_context = {
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
            for ls in line_segments],
        'knots': knots,
        'chords': chords,
        'rsft_generators': rsft_generators,
        'n_disks': len(disk_corners),
        'disk_corners': disk_corners,
        'x_labels': [{"label": x, "x": int(increment + increment * x + increment / 2)} for x in range(n_segments)],
        'y_labels': [{"label": y, "y": int(pad + increment * y + increment / 2)} for y in range(n_strands - 1)],
        'lch_graded_by': str(pd.lch_graded_by),
        'rsft_graded_by': str(pd.rsft_graded_by)
    }
    return template_context


@APP.route('/')
def home():
    data = request.args.get('link')
    if data is not None:
        pd = ll.PlatDiagram(**data)
    else:
        n_strands = request.args.get('n_strands')
        if n_strands is not None:
            n_strands = int(n_strands)
            crossings = request.args.get('crossings')
            if crossings is not None:
                crossings = [int(x) for x in crossings.split(",")]
            pd = ll.PlatDiagram(n_strands=n_strands, front_crossings=crossings)
        else:
            data = LINKS[DEFAULT_LINK_KEY].copy()
            data.pop('comment')
            pd = ll.PlatDiagram(**data)
    template_context = get_template_context(pd=pd)
    context = INDEX_TEMPLATE.render(**template_context)
    return make_response(context)


if __name__ == "__main__":
    APP.run()
