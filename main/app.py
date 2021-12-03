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

# Static resources
STATIC_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'static')
HTML_PATH = os.path.join(STATIC_PATH, 'index.html')
with open(HTML_PATH) as f:
    INDEX_TEMPLATE = Template(f.read())

# App resources
APP = Flask(__name__)


@APP.route('/')
def home():
    link = LINKS['m10_132_nv']
    del link['comment']
    pd = ll.PlatDiagram(**link)
    template_context = pd.get_svg_context()
    template_context['disk_corners'] = pd.disk_corners
    context = INDEX_TEMPLATE.render(**template_context)
    return make_response(context)


if __name__ == "__main__":
    APP.run()
