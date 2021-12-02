import utils


LOG = utils.get_logger(__name__)


class DiskSegment(object):
    def __init__(self, x, top_crossing_height=None, bottom_crossing_height=None,
                 left_crossing_height=None, right_crossing_height=None,
                 left_endpoints=None, right_endpoints=None):
        self.x = x
        self.top_crossing_height = top_crossing_height
        self.bottom_crossing_height = bottom_crossing_height
        self.left_crossing_height = left_crossing_height
        self.right_crossing_height = right_crossing_height
        self.left_endpoints = left_endpoints
        self.right_endpoints = right_endpoints

    def to_dict(self):
        return {
            "top_crossing_height": self.top_crossing_height,
            "bottom_crossing_height": self.bottom_crossing_height,
            "left_crossing_height": self.left_crossing_height,
            "right_crossing_height": self.right_crossing_height,
            "left_endpoints": self.left_endpoints,
            "right_endpoints": self.right_endpoints
        }


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
        if (vertex.x == 0) or (vertex.left_crossing_height is not None):
            return True
        return False

    def vertex_is_terminal(self, i):
        vertex = self.vertices[i]
        if (vertex.x == self.n_segments - 1) or (vertex.right_crossing_height is not None):
            return True
        return False

    def vertex_dead_end(self, i):
        for e in self.edges:
            if e[0] == i:
                return True
        return False

    def vertices_are_adjacent(self, i, j):
        v_left = self.vertices[i]
        v_right = self.vertices[j]
        if v_right.x != v_left.x + 1:
            return False
        return v_left.right_endpoints == v_right.left_endpoints

    def compute_edges(self):
        for i in range(len(self.vertices)):
            for j in range(len(self.vertices)):
                if self.vertices_are_adjacent(i, j):
                    self.edges.append([i, j])

    #

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
        # paths = [p for p in paths if self.vertex_is_terminal(p[-1])]
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

    def get_asymptotics(self):
        #TODO: This function appears to be broken
        asymptotics = []
        for d in self:
            if d.left_crossing_height is not None:
                asymptotics.append({
                    'x': d.x,
                    'y': d.left_crossing_height,
                    'z': 1
                })
            if d.bottom_crossing_height is not None:
                asymptotics.append({
                    'x': d.x,
                    'y': d.bottom_crossing_height,
                    'z': -1
                })
        for d in reversed(self):
            if d.right_crossing_height is not None:
                asymptotics.append({
                    'x': d.x,
                    'y': d.right_crossing_height,
                    'z': 1
                })
            if d.top_crossing_height is not None:
                asymptotics.append({
                    'x': d.x,
                    'y': d.top_crossing_height,
                    'z': -1
                })
        return asymptotics


class PlatSegment(object):
    def __init__(self, x, n_strands, crossing_height=None, left_close=False, right_close=False):
        self.x = x
        if n_strands % 2 != 0:
            raise RuntimeError("n_strands should be even")
        if n_strands < 1:
            raise RuntimeError("n_strands must be at least 2")
        self.n_strands = n_strands
        self.left_close = left_close
        self.right_close = right_close
        self.crossing_height = crossing_height

        self.left_knot_labels = [None for _ in range(n_strands)]
        self.right_knot_labels = [None for _ in range(n_strands)]

    def set_left_knot_label(self, i, j):
        self.left_knot_labels[i] = j

    def set_right_knot_label(self, i, j):
        self.left_knot_labels[i] = j

    def get_disk_segments(self):
        # enumerate all of the disk segments and return them as a set
        disk_segments = []
        # if it is a left close of a right close
        if self.left_close:
            for i in range(0, self.n_strands):
                if i % 2 == 0:
                    disk_segments.append(DiskSegment(x=self.x, right_endpoints=[i, i+1]))
            return disk_segments
        if self.right_close:
            for i in range(0, self.n_strands):
                if i % 2 == 0:
                    disk_segments.append(DiskSegment(x=self.x, left_endpoints=[i, i+1]))
            return disk_segments

        # the ones without crossings are
        for i in range(0, self.n_strands - 1):
            for j in range(i + 1, self.n_strands):
                if self.crossing_height not in [i, j]:
                    disk_segments.append(DiskSegment(
                        x=self.x, left_endpoints=[i, j], right_endpoints=[i, j]))
        # with an upward shift on the bottom
        for top_height in range(self.crossing_height + 2, self.n_strands):
            disk_segments.append(DiskSegment(
                x=self.x,
                left_endpoints=[self.crossing_height, top_height],
                right_endpoints=[self.crossing_height + 1, top_height]))
        # with an upward shift on the top
        for bottom_height in range(0, self.crossing_height):
            disk_segments.append(DiskSegment(
                x=self.x,
                left_endpoints=[bottom_height, self.crossing_height],
                right_endpoints=[bottom_height, self.crossing_height + 1]))
        # with a downward shift on the bottom
        for top_height in range(self.crossing_height + 2, self.n_strands):
            disk_segments.append(DiskSegment(
                x=self.x,
                left_endpoints=[self.crossing_height + 1, top_height],
                right_endpoints=[self.crossing_height, top_height]))
        # with a downward shift on the top
        for bottom_height in range(0, self.crossing_height):
            disk_segments.append(DiskSegment(
                x=self.x,
                left_endpoints=[bottom_height, self.crossing_height + 1],
                right_endpoints=[bottom_height, self.crossing_height]))
        # with a crossing on the top
        for bottom_height in range(0, self.crossing_height):
            disk_segments.append(DiskSegment(
                x=self.x,
                left_endpoints=[bottom_height, self.crossing_height],
                right_endpoints=[bottom_height, self.crossing_height],
                top_crossing_height=self.crossing_height))
        # with a crossing on the bottom
        for top_height in range(self.crossing_height + 1, self.crossing_height):
            disk_segments.append(DiskSegment(
                x=self.x,
                left_endpoints=[self.crossing_height + 1, top_height],
                right_endpoints=[self.crossing_height + 1, top_height],
                bottom_crossing_height=self.crossing_height))
        # with a crossing on the left
        disk_segments.append(DiskSegment(
            x=self.x,
            right_endpoints=[self.crossing_height, self.crossing_height + 1],
            left_crossing_height=self.crossing_height
        ))
        # with a crossing on the right
        disk_segments.append(DiskSegment(
            x=self.x,
            left_endpoints=[self.crossing_height, self.crossing_height + 1],
            right_crossing_height=self.crossing_height
        ))
        return disk_segments


class PlatDiagram(object):

    def __init__(self, n_strands, crossings=None):
        self.n_strands = n_strands

        # build self.plat_segments
        x = 0
        self.plat_segments = [PlatSegment(x=x, n_strands=self.n_strands, left_close=True)]
        x = 1
        # add internal segments with crossings in the style of Sivek's software
        if crossings is not None:
            for c in crossings:
                self.plat_segments.append(PlatSegment(x=x, n_strands=self.n_strands, crossing_height=c))
                x += 1
        # add crossings for right-pointing cusps
        for i in range(n_strands):
            if i % 2 == 0:
                self.plat_segments.append(PlatSegment(x=x, n_strands=self.n_strands, crossing_height=i))
                x += 1
        self.plat_segments.append(PlatSegment(x=x, n_strands=self.n_strands, right_close=True))

        self.crossings = self._get_crossings()
        self.disk_graph = self._get_disk_graph()
        self.disks = self.disk_graph.compute_disks()
        LOG.info(f"Disks in plat diagram: {len(self.disks)}")
        self.disk_asymptotics = self._get_disk_asymptotics()

    def _get_crossings(self):
        return [{"x": x, "y": self.plat_segments[x].crossing_height}
                for x in range(1, len(self.plat_segments) - 1)]

    def _get_disk_asymptotics(self):
        return [d.get_asymptotics() for d in self.disks]

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

    def _add_knot_labels_from_left(self, height, knot_index):
        self.plat_segments[0].set_left_knot_label(height, knot_index)
        self.plat_segments[0].set_right_knot_label(height, knot_index)
        current_height = 0
        for i in range(1, len(self.plat_segments)):
            self.plat_segments[i].set_left_knot_label(current_height, knot_index)
            if self.plat_segments[i].crossing_height == current_height:
                current_height += 1
            if self.plat_segments[i].crossing_height == current_height - 1:
                current_height -= 1
            self.plat_segments[i].set_right_knot_label(current_height, knot_index)

        self.plat_segments[-1].set_left_knot_label(current_height, knot_index)
