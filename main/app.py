import json
import os
from flask import (
    Flask,
    make_response,
    request,
    redirect,
    url_for)
from jinja2 import Template
import legendrian_links as ll

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
    pd = ll.PlatDiagram(4, [1, 1, 1])
    context = INDEX_TEMPLATE.render(
        plat_svg=pd.get_svg(),
        disk_corners=[[dc.to_dict() for dc in d] for d in pd.disk_corners])
    return make_response(context)


if __name__ == "__main__":
    APP.run()
