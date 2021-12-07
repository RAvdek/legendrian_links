import os
from jinja2 import Template
import utils


LOG = utils.get_logger(__name__)


class Crossing(object):

    def __init__(self, x, y):
        self.x = x
        self.y = y


class DiskCorner(object):
    CORNER_TYPES = {
        'l': '+',
        'r': '+',
        'u': '-',
        'd': '-'
    }

    def __init__(self, crossing, corner):
        self.crossing = crossing
        if corner not in self.CORNER_TYPES.keys():
            raise ValueError(f"Incoming type {corner} not in {self.CORNER_TYPES}")
        self.corner = corner

    def to_dict(self):
        return {
            'x': self.crossing.x,
            'y': self.crossing.y,
            'corner': self.corner,
            'pos_neg': self.CORNER_TYPES[self.corner]
        }


class DiskSegment(object):
    def __init__(self, x, disk_corner, left_endpoints=None, right_endpoints=None):
        self.x = x
        self.disk_corner = disk_corner
        self.corner = None
        if disk_corner is not None:
            self.corner = disk_corner.corner
        self.left_endpoints = left_endpoints
        self.right_endpoints = right_endpoints

    def get_endpoints(self):
        return [self.left_endpoints, self.right_endpoints]

    def to_dict(self):
        output = dict()
        if self.disk_corner is not None:
            output = self.disk_corner.to_dict()
        output['x'] = self.x
        output['left_endpoints'] = self.left_endpoints
        output['right_endpoints'] = self.right_endpoints
        return output


class DiskSegmentGraph(object):
    """A directed graph which encodes how DiskSegments with x coordinates are linked together"""

    def __init__(self, n_segments):
        self.vertices = list()
        self.edges = list()
        self.n_segments = n_segments

    def add_vertex(self, disk_segment):
        self.vertices.append(disk_segment)

    def vertex_is_initial(self, i):
        vertex = self.vertices[i]
        if (vertex.x == 0) or (vertex.corner == 'l'):
            return True
        return False

    def vertex_is_terminal(self, i):
        vertex = self.vertices[i]
        if (vertex.x == self.n_segments - 1) or (vertex.corner == 'r'):
            return True
        return False

    def vertex_dead_end(self, i):
        for e in self.edges:
            if e[0] == i:
                return True
        return False

    def vertices_are_consecutive(self, i, j):
        v_left = self.vertices[i]
        v_right = self.vertices[j]
        if v_right.x != v_left.x + 1:
            return False
        return v_left.right_endpoints == v_right.left_endpoints

    def compute_edges(self):
        for i in range(len(self.vertices)):
            for j in range(len(self.vertices)):
                if self.vertices_are_consecutive(i, j):
                    self.edges.append([i, j])

    def compute_paths_from_vertex(self, i):
        if self.vertex_is_terminal(i):
            return [[i]]
        outgoing_vertices = [e[1] for e in self.edges if e[0] == i]
        return [[i] + p for v in outgoing_vertices for p in self.compute_paths_from_vertex(v)]

    def compute_paths(self):
        paths = []
        for i in range(len(self.vertices)):
            if self.vertex_is_initial(i):
                paths += self.compute_paths_from_vertex(i)
        paths = [p for p in paths if self.vertex_is_terminal(p[-1])]
        LOG.info(f"DG paths: {paths}")
        return paths

    def path_to_disk(self, index_path):
        return Disk([self.vertices[i] for i in index_path])

    def compute_disks(self):
        paths = self.compute_paths()
        return [self.path_to_disk(p) for p in paths]


class Disk(list):
    """A disk is a list of disk segments such that the right endpoints of one segment
    agree with the left endpoints of the next. This is supposed to model an index 1
    J-disk in the Lagrangian projection determined by a Lagrangian resolution."""

    def get_endpoints(self):
        return [ds.get_endpoints() for ds in self]

    def get_disk_corners(self):
        disk_corners = []
        for ds in self:
            if ds.disk_corner is not None:
                if ds.disk_corner.corner in ['l', 'd']:
                    disk_corners.append(ds.disk_corner)
        for ds in reversed(self):
            if ds.disk_corner is not None:
                if ds.disk_corner.corner in ['u', 'r']:
                    disk_corners.append(ds.disk_corner)
        return disk_corners


class PlatSegment(object):
    def __init__(self, x, n_strands, crossing_y=None, left_close=False, right_close=False):
        self.x = x
        if n_strands % 2 != 0:
            raise RuntimeError("n_strands should be even")
        if n_strands < 1:
            raise RuntimeError("n_strands must be at least 2")
        self.n_strands = n_strands
        self.left_close = left_close
        self.right_close = right_close
        self.crossing_y = crossing_y

        self.left_knot_labels = [None]*n_strands
        self.right_knot_labels = [None]*n_strands

    def to_dict(self):
        return {
            "x": self.x,
            "n_strands": self.n_strands,
            "left_close": self.left_close,
            "right_close": self.right_close,
            "crossing_y": self.crossing_y
        }

    def get_disk_segments(self):
        # enumerate all of the disk segments and return them as a set
        disk_segments = []
        # if it is a left close of a right close
        if self.left_close:
            for i in range(0, self.n_strands):
                if i % 2 == 0:
                    disk_segments.append(DiskSegment(x=self.x, disk_corner=None, right_endpoints=[i, i + 1]))
            return disk_segments
        if self.right_close:
            for i in range(0, self.n_strands):
                if i % 2 == 0:
                    disk_segments.append(DiskSegment(x=self.x, disk_corner=None, left_endpoints=[i, i + 1]))
            return disk_segments

        # completely horizontal
        for i in range(0, self.n_strands - 1):
            for j in range(i + 1, self.n_strands):
                if self.crossing_y not in [i - 1, i, j - 1, j]:
                    disk_segments.append(DiskSegment(
                        x=self.x,
                        disk_corner=None,
                        left_endpoints=[i, j],
                        right_endpoints=[i, j]))
        # with an downward shift on the top
        for bottom_height in range(self.crossing_y + 2, self.n_strands):
            disk_segments.append(DiskSegment(
                x=self.x,
                disk_corner=None,
                left_endpoints=[self.crossing_y, bottom_height],
                right_endpoints=[self.crossing_y + 1, bottom_height]))
        # with an upward shift on the bottom
        for top_height in range(0, self.crossing_y):
            disk_segments.append(DiskSegment(
                x=self.x,
                disk_corner=None,
                left_endpoints=[top_height, self.crossing_y],
                right_endpoints=[top_height, self.crossing_y + 1]))
        # with a downward shift on the top
        for bottom_height in range(self.crossing_y + 1, self.n_strands):
            disk_segments.append(DiskSegment(
                x=self.x,
                disk_corner=None,
                left_endpoints=[self.crossing_y + 1, bottom_height],
                right_endpoints=[self.crossing_y, bottom_height]))
        # with a downward shift on the bottom
        for top_height in range(0, self.crossing_y):
            disk_segments.append(DiskSegment(
                x=self.x,
                disk_corner=None,
                left_endpoints=[top_height, self.crossing_y + 1],
                right_endpoints=[top_height, self.crossing_y]))
        # a crossing appears in all other cases
        crossing = Crossing(x=self.x, y=self.crossing_y)
        # with a crossing on the bottom
        for top_height in range(0, self.crossing_y):
            disk_segments.append(DiskSegment(
                x=self.x,
                disk_corner=DiskCorner(crossing=crossing, corner='d'),
                left_endpoints=[top_height, self.crossing_y],
                right_endpoints=[top_height, self.crossing_y]
            ))
        # with a crossing on the bottom
        for top_height in range(self.crossing_y + 1, self.n_strands):
            disk_segments.append(DiskSegment(
                x=self.x,
                disk_corner=DiskCorner(crossing=crossing, corner='u'),
                left_endpoints=[self.crossing_y + 1, top_height],
                right_endpoints=[self.crossing_y + 1, top_height],
            ))
        # with a crossing on the left
        disk_segments.append(DiskSegment(
            x=self.x,
            disk_corner=DiskCorner(crossing=crossing, corner='l'),
            right_endpoints=[self.crossing_y, self.crossing_y + 1],
        ))
        # with a crossing on the right
        disk_segments.append(DiskSegment(
            x=self.x,
            disk_corner=DiskCorner(crossing=crossing, corner='r'),
            left_endpoints=[self.crossing_y, self.crossing_y + 1],
        ))
        return disk_segments


class PlatDiagram(object):

    def __init__(self, n_strands, front_crossings=None):
        self.n_strands = n_strands
        self.front_crossings = front_crossings

        # build self.plat_segments
        x = 0
        self.plat_segments = [PlatSegment(x=x, n_strands=self.n_strands, left_close=True)]
        x = 1
        # add internal segments with crossings in the style of Sivek's software
        if front_crossings is not None:
            for c in front_crossings:
                self.plat_segments.append(PlatSegment(x=x, n_strands=self.n_strands, crossing_y=c))
                x += 1
        # add crossings for right-pointing cusps
        for i in range(n_strands):
            if i % 2 == 0:
                self.plat_segments.append(PlatSegment(x=x, n_strands=self.n_strands, crossing_y=i))
                x += 1
        self.plat_segments.append(PlatSegment(x=x, n_strands=self.n_strands, right_close=True))

        self.crossings = self._get_crossings()
        self.disk_graph = self._get_disk_graph()
        self.disks = self.disk_graph.compute_disks()
        LOG.info(f"Disks in plat diagram: {len(self.disks)}")
        self.disk_corners = self._get_disk_corners()

    def get_n_copy(self):
        #TODO
        pass

    def get_line_segments(self):
        n_segments = len(self.plat_segments)
        lines = []

        for ps in self.plat_segments:
            if ps.left_close:
                for i in range(self.n_strands):
                    if i % 2 == 0:
                        lines.append([[0, i + .5], [1, i]])
                    else:
                        lines.append([[0, i - .5], [1, i]])
            elif ps.right_close:
                for i in range(self.n_strands):
                    if i % 2 == 0:
                        lines.append([[n_segments - 1, i], [n_segments, i + .5]])
                    else:
                        lines.append([[n_segments - 1, i], [n_segments, i - .5]])
            else:
                x = ps.x
                crossing_y = ps.crossing_y
                for s in range(self.n_strands):
                    if s == crossing_y:
                        lines.append([[x, s], [x + 1, s + 1]])
                    elif s == crossing_y + 1:
                        lines.append([[x, s], [x + 1, s - 1]])
                    else:
                        lines.append([[x, s],[x + 1, s]])
        return lines


    def _get_crossings(self):
        return [Crossing(x=x, y=self.plat_segments[x].crossing_y)
                for x in range(1, len(self.plat_segments) - 1)]

    def _get_disk_corners(self):
        return [d.get_disk_corners() for d in self.disks]

    def _get_disk_graph(self):
        disk_graph = DiskSegmentGraph(n_segments=len(self.plat_segments))
        for x in range(len(self.plat_segments)):
            segment = self.plat_segments[x]
            disk_segments = segment.get_disk_segments()
            for d in disk_segments:
                disk_graph.add_vertex(d)
        disk_graph.compute_edges()
        LOG.info(f"{len(disk_graph.vertices)} disk_graph vertices")
        LOG.info(f"{len(disk_graph.edges)} disk_graph edges")
        return disk_graph

    def _assign_knot_labels_xy(self, knot_label, x=0, y=0):
        #TODO: Update method to include orientations
        left_close = self.plat_segments[x].left_close
        right_close = self.plat_segments[x].right_close
        crossing_height = self.plat_segments[x].crossing_y

        if left_close or right_close:
            self.plat_segments[x].left_knot_labels[y] = knot_label
            if y % 2 == 0:
                self.plat_segments[x].left_knot_labels[y + 1] = knot_label
                self.plat_segments[x].right_knot_labels[y + 1] = knot_label
            else:
                self.plat_segments[x].left_knot_labels[y - 1] = knot_label
                self.plat_segments[x].right_knot_labels[y - 1] = knot_label
        else:
            self.plat_segments[x].left_knot_labels[y] = knot_label
            if crossing_height == y:
                self.plat_segments[x].right_knot_labels[y + 1] = knot_label
            elif crossing_height == y-1:
                self.plat_segments[x].right_knot_labels[y - 1] = knot_label
            else:
                self.plat_segments[x].right_knot_labels[y] = knot_label
