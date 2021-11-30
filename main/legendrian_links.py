import utils


LOG = utils.get_logger(__name__)


class DiskSegment(object):
    def __init__(self, top_crossing_height=None, bottom_crossing_height=None,
                 left_crossing_height=None, right_crossing_height=None,
                 left_endpoints=None, right_endpoints=None):
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


class DiskGraph(object):
    """A directed graph which encodes how DiskSegments with x coordinates are linked together"""

    def __init__(self, n_segments):
        self.vertices = list()
        self.edges = list()
        self.n_segments = n_segments

    def add_vertex(self, x, disk_segment):
        self.vertices.append({
            "x": x,
            "disk_segment": disk_segment
        })

    def vertex_is_initial(self, i):
        vertex = self.vertices[i]
        if (vertex["x"] == 0) or (vertex["disk_segment"].left_crossing_height is not None):
            return True
        return False

    def vertex_is_terminal(self, i):
        vertex = self.vertices[i]
        if (vertex["x"] == self.n_segments) or (vertex["disk_segment"].right_crossing_height is not None):
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
        if v_right["x"] != v_left["x"] + 1:
            return False
        return v_left["disk_segment"].right_endpoints == v_right["disk_segment"].left_endpoints

    def compute_edges(self):
        for i in range(len(self.vertices)):
            for j in range(len(self.vertices)):
                if self.vertices_are_adjacent(i, j):
                    self.edges.append([i, j])

    def compute_paths_from_vertex(self, i):
        if self.vertex_is_terminal(i):
            return [i]
        outgoing_vertices = [e[1] for e in self.edges if e[0] == i]
        return [[i] + self.compute_paths_from_vertex(v) for v in outgoing_vertices]

    def compute_paths(self):
        paths = []
        for i in range(len(self.vertices)):
            if self.vertex_is_initial(i):
                paths += self.compute_paths_from_vertex(i)
        paths = [p for p in paths if self.vertex_is_terminal(p[-1])]
        return paths

    def path_to_disk(self, index_path):
        return [self.vertices[i]["disk_segment"] for i in index_path]

    def compute_disks(self):
        paths = self.compute_paths()
        return [self.path_to_disk(p) for p in paths]


class Disk(list):
    """A disk is a list of disk segments such that the right endpoints of one segment
    agree with the left endpoints of the next. This is supposed to model an index 1
    J-disk in the Lagrangian projection determined by a Lagrangian resolution."""


class PlatSegment(object):
    def __init__(self, n_strands, crossing_height=None, left_close=False, right_close=False):
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
                    disk_segments.append(DiskSegment(
                        right_endpoints=[i, i+1]
                    ))
            return disk_segments
        if self.right_close:
            for i in range(0, self.n_strands):
                if i % 2 == 0:
                    disk_segments.append(DiskSegment(
                        left_endpoints=[i, i+1]
                    ))
            return disk_segments

        # the ones without crossings are
        for i in range(0, self.n_strands - 1):
            for j in range(i + 1, self.n_strands):
                if self.crossing_height not in [i, j]:
                    disk_segments.append(DiskSegment(
                        left_endpoints=[i, j], right_endpoints=[i, j]))
        # with an upward shift on the bottom
        for top_height in range(self.crossing_height + 2, self.n_strands):
            disk_segments.append(DiskSegment(
                left_endpoints=[self.crossing_height, top_height],
                right_endpoints=[self.crossing_height + 1, top_height]))
        # with an upward shift on the top
        for bottom_height in range(0, self.crossing_height):
            disk_segments.append(DiskSegment(
                left_endpoints=[bottom_height, self.crossing_height],
                right_endpoints=[bottom_height, self.crossing_height + 1]))
        # with a downward shift on the bottom
        for top_height in range(self.crossing_height + 2, self.n_strands):
            disk_segments.append(DiskSegment(
                left_endpoints=[self.crossing_height + 1, top_height],
                right_endpoints=[self.crossing_height, top_height]))
        # with a downward shift on the top
        for bottom_height in range(0, self.crossing_height):
            disk_segments.append(DiskSegment(
                left_endpoints=[bottom_height, self.crossing_height + 1],
                right_endpoints=[bottom_height, self.crossing_height]))
        # with a crossing on the top
        for bottom_height in range(0, self.crossing_height):
            disk_segments.append(DiskSegment(
                left_endpoints=[bottom_height, self.crossing_height],
                right_endpoints=[bottom_height, self.crossing_height],
                top_crossing_height=self.crossing_height))
        # with a crossing on the bottom
        for top_height in range(self.crossing_height + 1, self.crossing_height):
            disk_segments.append(DiskSegment(
                left_endpoints=[self.crossing_height + 1, top_height],
                right_endpoints=[self.crossing_height + 1, top_height],
                bottom_crossing_height=self.crossing_height))
        # with a crossing on the left
        disk_segments.append(DiskSegment(
            right_endpoints=[self.crossing_height, self.crossing_height + 1],
            left_crossing_height=self.crossing_height
        ))
        # with a crossing on the right
        disk_segments.append(DiskSegment(
            left_endpoints=[self.crossing_height, self.crossing_height + 1],
            right_crossing_height=self.crossing_height
        ))
        return disk_segments


class PlatDiagram(object):

    def __init__(self, n_strands, crossings=None):
        self.n_strands = n_strands
        self.segments = [PlatSegment(n_strands=self.n_strands, left_close=True)]
        # add internal segments with crossings in the style of Sivek's software
        if crossings is not None:
            self.segments += [
                PlatSegment(n_strands=self.n_strands, crossing_height=x) for x in crossings]
        # add crossings for right-pointing cusps
        for i in range(n_strands):
            if i % 2 == 0:
                self.segments.append(PlatSegment(n_strands=self.n_strands, crossing_height=i))
        self.segments.append(PlatSegment(n_strands=self.n_strands, right_close=True))

    def del_inner_segment(self, i):
        self.segments.pop(i)

    def add_inner_segment(self, i, crossing=None):
        self.segments.insert(i, PlatSegment(n_strands=self.n_strands, crossing_height=crossing))

    def get_disk_graph(self):
        disk_graph = DiskGraph(n_segments=len(self.segments))
        for x in range(len(self.segments)):
            segment = self.segments[x]
            disk_segments = segment.get_disk_segments()
            for d in disk_segments:
                disk_graph.add_vertex(x, d)
        disk_graph.compute_edges()
        return disk_graph

    def get_disks(self):
        disk_graph = self.get_disk_graph()
        return disk_graph.compute_disks()

    def _add_knot_labels_from_left(self, height, knot_index):
        self.segments[0].set_left_knot_label(height, knot_index)
        self.segments[0].set_right_knot_label(height, knot_index)
        current_height = 0
        for i in range(1, len(self.segments)):
            self.segments[i].set_left_knot_label(current_height, knot_index)
            if self.segments[i].crossing_height == current_height:
                current_height += 1
            if self.segments[i].crossing_height == current_height - 1:
                current_height -= 1
            self.segments[i].set_right_knot_label(current_height, knot_index)

        self.segments[-1].set_left_knot_label(current_height, knot_index)
