import utils


LOG = utils.get_logger(__name__)


class Crossing(object):
    """Object which carries crossing information.
    Should be able to see the line segments at the positive and negative ends of a chord.
    From this we should be able to implement signs for holomorphic disks."""


    def __init__(self, top_line_segment, bottom_line_segment):
        self.top_line_segment = top_line_segment
        self.bottom_line_segment = bottom_line_segment


class DiskCorner(object):
    """Convex corner of a disk in the Lagrangian projection. Corner type tells the direction it points."""
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
            'x': self.crossing.top_line_segment.x_left,
            'y': self.crossing.top_line_segment.y_left,
            'from_knot': self.crossing.bottom_line_segment.knot_label,
            'to_knot': self.crossing.top_line_segment.knot_label,
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

    def __init__(self, line_segments, left_close=False, right_close=False):
        self.line_segments = line_segments
        self.n_strands = len(line_segments)
        self.left_close = left_close
        self.right_close = right_close
        self._set_x()
        self._set_crossing_data()

    def _set_x(self):
        x_left_values = list(set([ls.x_left for ls in self.line_segments]))
        n_x_left_values = len(x_left_values)
        if n_x_left_values > 1:
            raise ValueError(f"Found {n_x_left_values} different left x values for line segments in PlatSegment")
        self.x = x_left_values[0]

    def _set_crossing_data(self):
        if self.left_close or self.right_close:
            self.crossing_y = None
            self.top_line_segment_at_crossing = None
            self.bottom_line_segment_at_crossing = None
        else:
            southeast_segments = [ls for ls in self.line_segments if ls.y_left < ls.y_right]
            n_southeast_segments = len(southeast_segments)
            northeast_segments = [ls for ls in self.line_segments if ls.y_left > ls.y_right]
            n_northeast_segments = len(northeast_segments)
            if n_southeast_segments != 1 or n_northeast_segments != 1:
                raise ValueError(f"Trying to initialize PlatSegment at x={self.x} with "
                                 f"{n_southeast_segments} NE segments and "
                                 f"{n_northeast_segments} NE segments. "
                                 f"Line segments: {[ls.to_array() for ls in self.line_segments]}")
            self.crossing_y = southeast_segments[0].y_left
            self.top_line_segment_at_crossing = southeast_segments[0]
            self.bottom_line_segment_at_crossing = northeast_segments[0]

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
        crossing = Crossing(
            top_line_segment=self.top_line_segment_at_crossing,
            bottom_line_segment=self.bottom_line_segment_at_crossing)
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


class LineSegment(object):
    ORIENTATIONS = ['l', 'r']

    def __init__(self, x_left, y_left, x_right, y_right):
        self.x_left = x_left
        self.y_left = y_left
        self.x_right = x_right
        self.y_right = y_right
        self.orientation = None
        self.knot_label = None

    def to_array(self):
        return [[self.x_left, self.y_left], [self.x_right, self.y_right]]

    def set_orientation(self, o):
        if o not in self.ORIENTATIONS:
            raise ValueError(f'LineSegment orientation must be in {self.ORIENTATIONS}')
        self.orientation = o

    def set_knot_label(self, kl):
        self.knot_label = kl


class PlatDiagram(object):

    def __init__(self, n_strands, front_crossings=None):
        self.n_strands = n_strands
        self.front_crossings = front_crossings

        self._set_line_segments()
        self._label_line_segments()
        self._set_plat_segments()
        self._set_disk_graph()
        self._set_disks()
        LOG.info(f"Disks in plat diagram: {len(self.disks)}")
        self._set_disk_corners()


    def get_line_segment_arrays(self):
        return [ls.to_array() + [ls.knot_label, ls.orientation] for ls in self.line_segments]


    def _set_line_segments(self):
        lines = []

        # make the left closing part
        x = 0
        for i in range(self.n_strands):
            if i % 2 == 0:
                lines.append(LineSegment(x_left=0, y_left=i + .5, x_right=1, y_right=i))
            else:
                lines.append(LineSegment(x_left=0, y_left=i - .5, x_right=1, y_right=i))

        if self.front_crossings is not None:
            for fcy in self.front_crossings:
                x += 1
                for s in range(self.n_strands):
                    if s == fcy:
                        lines.append(LineSegment(x_left=x, y_left=s, x_right=x + 1, y_right=s + 1))
                    elif s == fcy + 1:
                        lines.append(LineSegment(x_left=x, y_left=s, x_right=x + 1, y_right=s - 1))
                    else:
                        lines.append(LineSegment(x_left=x, y_left=s, x_right=x + 1, y_right=s))

        # make crossings for each of the right-pointing cusps
        for i in [_ for _ in range(self.n_strands) if _ % 2 == 0]:
            x += 1
            for j in range(self.n_strands):
                if j == i:
                    lines.append(LineSegment(x_left=x, y_left=j, x_right=x + 1, y_right=j + 1))
                elif j == i + 1:
                    lines.append(LineSegment(x_left=x, y_left=j, x_right=x + 1, y_right=j - 1))
                else:
                    lines.append(LineSegment(x_left=x, y_left=j, x_right=x + 1, y_right=j))

        # make the right closing part
        x += 1
        for i in range(self.n_strands):
            if i % 2 == 0:
                lines.append(LineSegment(x_left=x, y_left=i, x_right=x + 1, y_right=i + .5))
            else:
                lines.append(LineSegment(x_left=x, y_left=i, x_right=x + 1, y_right=i - .5))

        self.line_segments = lines


    def _label_line_segments(self):
        knot_label = 0
        # The while loops are in danger of being infinite. Cut at an n to debug.
        n = 0
        while True:
            unlabeled_line_segments = [ls for ls in self.line_segments if ls.knot_label is None]
            n_unlabeled_line_segments = len(unlabeled_line_segments)
            LOG.info(f"{n_unlabeled_line_segments} line segments unlabeled")
            if n_unlabeled_line_segments == 0:
                break
            initial_ls = unlabeled_line_segments[0]
            initial_ls.set_knot_label(knot_label)
            initial_ls.set_orientation('r')
            ls = self._label_next_line_segment(initial_ls)
            while ls != initial_ls:
                LOG.info(ls.to_array())
                ls = self._label_next_line_segment(ls)
                n += 1
            knot_label += 1
        self.n_components = len(set([ls.knot_label for ls in self.line_segments]))

    def _set_plat_segments(self):
        max_left_x = max([ls.x_left for ls in self.line_segments])
        plat_segments = []
        for x in range(0, max_left_x + 1):
            left_x_segments = [ls for ls in self.line_segments if ls.x_left == x]
            left_close = x == 0
            right_close = x == max_left_x
            ps = PlatSegment(left_x_segments, left_close=left_close, right_close=right_close)
            plat_segments.append(ps)
        self.plat_segments = plat_segments

    def _set_disk_graph(self):
        disk_graph = DiskSegmentGraph(n_segments=len(self.plat_segments))
        for x in range(len(self.plat_segments)):
            segment = self.plat_segments[x]
            disk_segments = segment.get_disk_segments()
            for d in disk_segments:
                disk_graph.add_vertex(d)
        disk_graph.compute_edges()
        LOG.info(f"{len(disk_graph.vertices)} disk_graph vertices")
        LOG.info(f"{len(disk_graph.edges)} disk_graph edges")
        self.disk_graph = disk_graph

    def _set_disks(self):
        self.disks = self.disk_graph.compute_disks()


    def _set_disk_corners(self):
        self.disk_corners = [d.get_disk_corners() for d in self.disks]


    def get_line_segment_by_right_xy(self, x, y):
        search_results = [ls for ls in self.line_segments if ls.x_right == x and ls.y_right == y]
        n_results = len(search_results)
        if n_results != 1:
            raise RuntimeError(
                f'Found {n_results} in search for line segment by right_xy, rather than 1 at x,y={x},{y}. '
                f'Line segments available: {[ls.to_array() for ls in self.line_segments]}')
        return search_results[0]


    def get_line_segment_by_left_xy(self, x, y):
        search_results = [ls for ls in self.line_segments if ls.x_left == x and ls.y_left == y]
        n_results = len(search_results)
        if n_results != 1:
            raise RuntimeError(
                f'Found {n_results} in search for line segment by left_xy, rather than 1 at x,y={x},{y}. '
                f'Line segments available: {[ls.to_array() for ls in self.line_segments]}')
        return search_results[0]


    def _label_next_line_segment(self, ls):
        if ls not in self.line_segments:
            raise ValueError('line_segment does not belong to PlatDiagram instance')
        if ls.orientation is None:
            raise ValueError('Cannot find next line segment if it is not orientated')
        if ls.knot_label is None:
            raise ValueError('Cannot find knot_label for line_segment')

        max_x_right = max([ls.x_right for ls in self.line_segments])
        if ls.orientation == 'r':
            if ls.x_right < max_x_right:
                next_ls = self.get_line_segment_by_left_xy(x=ls.x_right, y=ls.y_right)
                next_ls.set_orientation('r')
            else:
                if ls.y_right < ls.y_left:
                    next_ls = self.get_line_segment_by_left_xy(x=ls.x_left, y=ls.y_left - 1)
                else:
                    next_ls = self.get_line_segment_by_left_xy(x=ls.x_left, y=ls.y_left + 1)
                next_ls.set_orientation('l')
        else:
            if ls.x_left > 0:
                next_ls = self.get_line_segment_by_right_xy(x=ls.x_left, y=ls.y_left)
                next_ls.set_orientation('l')
            else:
                if ls.y_right < ls.y_left:
                    next_ls = self.get_line_segment_by_right_xy(x=ls.x_right, y=ls.y_right + 1)
                else:
                    next_ls = self.get_line_segment_by_right_xy(x=ls.x_right, y=ls.y_right - 1)
                next_ls.set_orientation('r')

        next_ls.set_knot_label(ls.knot_label)
        return next_ls
