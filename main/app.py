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


def get_template_context(plat_diagram, increment=50, pad=10):
    n_strands = plat_diagram.n_strands
    n_segments = len(plat_diagram.plat_segments)
    height = increment * plat_diagram.n_strands + pad
    width = (increment * n_segments) + pad + increment
    line_segments = plat_diagram.get_line_segment_arrays()
    disk_corners = plat_diagram.disk_corners
    template_context = {
        "pad": pad,
        "increment": increment,
        "height": height,
        "width": width,
        "lines": [
            [
                [increment + pad + int(increment * l[0][0]), pad + int(increment * l[0][1])],
                [increment + pad + int(increment * l[1][0]), pad + int(increment * l[1][1])],
                KNOT_COLORS[l[2]],
            ]
            for l in line_segments],
        'n_disks': len(disk_corners),
        'disk_corners': disk_corners,
        'x_labels': [{"label": x, "x": int(increment + increment * x + increment / 2)} for x in range(n_segments)],
        'y_labels': [{"label": y, "y": int(pad + increment * y + increment / 2)} for y in range(n_strands - 1)]
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
            data = LINKS[DEFAULT_LINK_KEY]
            data.pop('comment')
            pd = ll.PlatDiagram(**data)
    template_context = get_template_context(plat_diagram=pd)
    context = INDEX_TEMPLATE.render(**template_context)
    return make_response(context)


if __name__ == "__main__":
    APP.run()
