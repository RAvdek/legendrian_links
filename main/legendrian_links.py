from collections import Counter
import utils


LOG = utils.get_logger(__name__)


class Chord(object):
    """Object which carries crossing information.
    Should be able to see the line segments at the positive and negative ends of a chord.
    From this we should be able to implement signs for holomorphic disks."""


    def __init__(self, top_line_segment, bottom_line_segment):
        self.top_line_segment = top_line_segment
        self.bottom_line_segment = bottom_line_segment
        self._set_xy()
        self._set_sign()
        self._set_grading()

    def to_string(self):
        return f"[{self.x}]"

    def is_composable_with(self, chord):
        return self.top_line_segment.knot_label == chord.bottom_line_segment.knot_label

    def _set_xy(self):
        self.x = self.top_line_segment.x_left
        self.y = self.top_line_segment.y_left

    def _set_sign(self):
        if self.top_line_segment.orientation is None or self.bottom_line_segment.orientation is None:
            raise ValueError(f"Trying to initialize crossing with unoriented line segments. "
                             f"Top {self.top_line_segment.to_array()}, "
                             f"bottom {self.bottom_line_segment.to_array()}")
        self.sign = 1 if self.top_line_segment.orientation == self.bottom_line_segment.orientation else -1

    def _set_grading(self):
        self.grading = 1 if self.sign == -1 else 0


class WordOfChords(object):

    def __init__(self, word):
        for chord in word:
            if not isinstance(chord, Chord):
                raise ValueError("Trying to create WordOfChord object not from a list of chords")
        self.word = word
        self._set_grading()
        self._set_x_array()

    def to_string(self):
        return "[" + ",".join([str(chord.x) for chord in self.word]) + "]"

    def _set_grading(self):
        self.grading = (-1 + len(self.word) + len([chord for chord in self.word if chord.sign == -1])) % 2

    def _set_x_array(self):
        self.x_array = [chord.bottom_line_segment.x_left for chord in self.word]


class DiskCorner(object):
    """Convex corner of a disk in the Lagrangian projection. Corner type tells the direction it points."""
    CORNER_TYPES = {
        'l': '+',
        'r': '+',
        'u': '-',
        'd': '-'
    }

    def __init__(self, chord, corner):
        self.chord = chord
        if corner not in self.CORNER_TYPES.keys():
            raise ValueError(f"Incoming type {corner} not in {self.CORNER_TYPES}")
        self.corner = corner
        self.pos_neg = self.CORNER_TYPES[self.corner]

    def to_dict(self):
        return {
            'x': self.chord.top_line_segment.x_left,
            'y': self.chord.top_line_segment.y_left,
            'from_knot': self.chord.bottom_line_segment.knot_label,
            'to_knot': self.chord.top_line_segment.knot_label,
            'corner': self.corner,
            'pos_neg': self.pos_neg
        }


class DiskSegment(object):
    # TODO: rather than having left and right endpoints, it might be easier to have left_top,
    #  left_bottom, right_top, right_bottom line segments. Current implementation is a hack.

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
        # TODO: Divide and conquer is suboptimal. Can instead use dynamic which is faster.
        # This may actually be important for very large numbers of disks.
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


class Disk(object):
    """A disk is a list of disk segments such that the right endpoints of one segment
    agree with the left endpoints of the next. This is supposed to model an index 1
    J-disk in the Lagrangian projection determined by a Lagrangian resolution."""

    def __init__(self, disk_segments):
        self.disk_segments = disk_segments
        self._set_disk_corners()

    def get_endpoints(self):
        return [ds.get_endpoints() for ds in self.disk_segments]

    def _set_disk_corners(self):
        disk_corners = []
        disk_segments = self.disk_segments
        for ds in disk_segments:
            if ds.disk_corner is not None:
                if ds.disk_corner.corner in ['l', 'd']:
                    disk_corners.append(ds.disk_corner)
        for ds in reversed(disk_segments):
            if ds.disk_corner is not None:
                if ds.disk_corner.corner in ['u', 'r']:
                    disk_corners.append(ds.disk_corner)
        self.disk_corners = disk_corners
        self.pos_corners = [c for c in self.disk_corners if c.pos_neg == '+']
        self.neg_corners = [c for c in self.disk_corners if c.pos_neg == '-']

    def is_lch(self):
        return len(self.pos_corners) == 1

    def is_rsft(self):
        raise NotImplementedError()


class PlatSegment(object):

    def __init__(self, line_segments, chord=None, left_close=False, right_close=False):
        self.line_segments = line_segments
        self.chord = chord
        self.n_strands = len(line_segments)
        self.left_close = left_close
        self.right_close = right_close
        self._set_xy()

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
        chord = self.chord
        # with a crossing on the bottom
        for top_height in range(0, self.crossing_y):
            disk_segments.append(DiskSegment(
                x=self.x,
                disk_corner=DiskCorner(chord=chord, corner='d'),
                left_endpoints=[top_height, self.crossing_y],
                right_endpoints=[top_height, self.crossing_y]
            ))
        # with a crossing on the bottom
        for top_height in range(self.crossing_y + 1, self.n_strands):
            disk_segments.append(DiskSegment(
                x=self.x,
                disk_corner=DiskCorner(chord=chord, corner='u'),
                left_endpoints=[self.crossing_y + 1, top_height],
                right_endpoints=[self.crossing_y + 1, top_height],
            ))
        # with a crossing on the left
        disk_segments.append(DiskSegment(
            x=self.x,
            disk_corner=DiskCorner(chord=chord, corner='l'),
            right_endpoints=[self.crossing_y, self.crossing_y + 1],
        ))
        # with a crossing on the right
        disk_segments.append(DiskSegment(
            x=self.x,
            disk_corner=DiskCorner(chord=chord, corner='r'),
            left_endpoints=[self.crossing_y, self.crossing_y + 1],
        ))
        return disk_segments

    def _set_xy(self):
        x_left_values = list(set([ls.x_left for ls in self.line_segments]))
        n_x_left_values = len(x_left_values)
        if n_x_left_values > 1:
            raise ValueError(f"Found {n_x_left_values} different left x values for line segments in PlatSegment")
        self.x = x_left_values[0]
        self.crossing_y = None
        if self.chord is not None:
            self.crossing_y = self.chord.top_line_segment.y_left


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

    def __init__(self, n_strands, front_crossings=[]):
        self.n_strands = n_strands
        self.front_crossings = front_crossings

        self._set_line_segments()
        self.max_x_left = max([ls.x_left for ls in self.line_segments])
        self._label_line_segments()
        # wait to set crossings until the line segments have been labeled
        self._set_chords()
        self._set_knots()
        self._set_plat_segments()
        self._set_disk_graph()
        self._set_disks()
        LOG.info(f"Disks in plat diagram: {len(self.disks)}")
        self._set_disk_corners()
        self._set_lch_generators()
        self._set_lch_del()
        self._set_rsft_generators()

    def get_line_segment_arrays(self):
        return [ls.to_array() + [ls.knot_label, ls.orientation] for ls in self.line_segments]

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

    def get_next_line_segment(self, line_segment):
        if line_segment.orientation == 'r':
            if line_segment.x_left == self.max_x_left:
                if line_segment.y_left > line_segment.y_right:
                    return self.get_line_segment_by_left_xy(x=line_segment.x_left, y=line_segment.y_left - 1)
                else:
                    return self.get_line_segment_by_left_xy(x=line_segment.x_left, y=line_segment.y_left + 1)
            else:
                return self.get_line_segment_by_left_xy(x=line_segment.right_x, y=line_segment.right_y)
        else:
            if line_segment.x_left == 0:
                if line_segment.y_left > line_segment.y_right:
                    return self.get_line_segment_by_right_xy(x=line_segment.x_right, y=line_segment.y_right + 1)
                else:
                    return self.get_line_segment_by_right_xy(x=line_segment.x_right, y=line_segment.y_right - 1)
            else:
                return self.get_line_segment_by_right_xy(x=line_segment.left_x, y=line_segment.left_y)

    def line_segment_is_incoming_to_left_up_cusp(self, line_segment):
        if not (line_segment.x_left == 0 and line_segment.orientation == 'l'):
            return False
        next_ls = self.get_next_line_segment(line_segment)
        return next_ls.y_right < line_segment.y_right

    def line_segment_is_incoming_to_left_down_cusp(self, line_segment):
        if not (line_segment.x_left == 0 and line_segment.orientation == 'l'):
            return False
        next_ls = self.get_next_line_segment(line_segment)
        return next_ls.y_right > line_segment.y_right

    def line_segment_is_incoming_to_right_up_cusp(self, line_segment):
        if not (line_segment.x_left == self.max_x_left and line_segment.orientation == 'r'):
            return False
        next_ls = self.get_next_line_segment(line_segment)
        return next_ls.y_left < line_segment.y_left

    def line_segment_is_incoming_to_right_down_cusp(self, line_segment):
        if not (line_segment.x_left == self.max_x_left and line_segment.orientation == 'r'):
            return False
        next_ls = self.get_next_line_segment(line_segment)
        return next_ls.y_left > line_segment.y_left

    def _set_line_segments(self):
        lines = []
        lag_crossings_from_right_cusps = [i for i in range(self.n_strands) if i % 2 == 0]
        lag_crossings = self.front_crossings + lag_crossings_from_right_cusps

        # make the left closing part
        x = 0
        for i in range(self.n_strands):
            if i % 2 == 0:
                lines.append(LineSegment(x_left=0, y_left=i + .5, x_right=1, y_right=i))
            else:
                lines.append(LineSegment(x_left=0, y_left=i - .5, x_right=1, y_right=i))
        # make the segments at x values where there are crossings
        if lag_crossings is not None:
            for fcy in lag_crossings:
                x += 1
                for s in range(self.n_strands):
                    if s == fcy:
                        lines.append(LineSegment(x_left=x, y_left=s, x_right=x + 1, y_right=s + 1))
                    elif s == fcy + 1:
                        lines.append(LineSegment(x_left=x, y_left=s, x_right=x + 1, y_right=s - 1))
                    else:
                        lines.append(LineSegment(x_left=x, y_left=s, x_right=x + 1, y_right=s))
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

    def _label_next_line_segment(self, ls):
        if ls not in self.line_segments:
            raise ValueError('line_segment does not belong to PlatDiagram instance')
        if ls.orientation is None:
            raise ValueError('Cannot find next line segment if it is not orientated')
        if ls.knot_label is None:
            raise ValueError('Cannot find knot_label for line_segment')

        max_x_right = self.max_x_left + 1
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

    def _set_chords(self):
        chords = []
        for x in range(0, self.max_x_left):
            x_segments = [ls for ls in self.line_segments if ls.x_left == x]
            for y in range(self.n_strands - 1):
                if x == 0:
                    top_left_ls = [ls for ls in x_segments if ls.y_right == y][0]
                    bottom_left_ls = [ls for ls in x_segments if ls.y_right == y + 1][0]
                else:
                    top_left_ls = [ls for ls in x_segments if ls.y_left == y][0]
                    bottom_left_ls = [ls for ls in x_segments if ls.y_left == y + 1][0]
                are_crossing = (top_left_ls.y_right == y + 1) and (bottom_left_ls.y_right == y)
                if are_crossing:
                    chords.append(Chord(top_line_segment=top_left_ls, bottom_line_segment=bottom_left_ls))
        self.chords = chords

    def _set_knots(self):
        knot_labels = sorted(list(set([ls.knot_label for ls in self.line_segments])))
        knots = list()
        for kl in knot_labels:
            knot_line_segments = [ls for ls in self.line_segments if ls.knot_label == kl]
            # compute tb as writhe in the Lagrangian projection
            self_chords = [c for c in self.chords
                           if c.top_line_segment.knot_label == kl and c.bottom_line_segment.knot_label == kl]
            tb = sum([c.sign for c in self_chords])
            # compute rot from the left and right cusps
            rot = 0
            for ls in knot_line_segments:
                if self.line_segment_is_incoming_to_left_down_cusp(ls):
                    rot += 1
                if self.line_segment_is_incoming_to_right_up_cusp(ls):
                    rot += 1
                if self.line_segment_is_incoming_to_left_up_cusp(ls):
                    rot -= 1
                if self.line_segment_is_incoming_to_right_down_cusp(ls):
                    rot -= 1
            rot = rot // 2

            knots.append({
                'line_segments': knot_line_segments,
                'tb': tb,
                'rot': rot
            })

        self.knots = knots

    def _set_plat_segments(self):
        max_left_x = max([ls.x_left for ls in self.line_segments])
        plat_segments = []
        for x in range(0, max_left_x + 1):
            left_x_segments = [ls for ls in self.line_segments if ls.x_left == x]
            left_close = x == 0
            right_close = x == max_left_x
            chords = [c for c in self.chords
                      if c.bottom_line_segment in left_x_segments or c.top_line_segment in left_x_segments]
            if len(chords) > 1:
                raise ValueError(f'More than 1 chord at a plat segment')
            chord = chords[0] if len(chords) == 1 else None
            ps = PlatSegment(left_x_segments, chord = chord, left_close=left_close, right_close=right_close)
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
        self.disk_corners = [d.disk_corners for d in self.disks]

    def _set_lch_generators(self):
        self.lch_generators = self.chords

    def _set_lch_del(self):
        lch_disks = [d for d in self.disks if d.is_lch()]
        lch_del = {g: list() for g in self.lch_generators}
        for d in lch_disks:
            disk_corners = d.disk_corners.copy()
            # cyclically rotate disk until positive corner is in the last (-1) position
            while disk_corners[-1].pos_neg == '-':
                last_corner = disk_corners.pop(-1)
                disk_corners = [last_corner] + disk_corners
            pos_chord = disk_corners.pop().chord
            neg_chords = WordOfChords(word=[dc.chord for dc in disk_corners])
            lch_del[pos_chord].append(neg_chords)
        self.lch_del = lch_del

    def _set_rsft_generators(self):
        # length 1 composable words
        composable_admissible_words = [[[chord] for chord in self.chords]]
        # dynamic algo to find words of length n from words of length n-1
        for _ in range(2, len(self.knots)+1):
            largest_length_words = composable_admissible_words[-1]
            new_words = [
                word + [chord] for word in largest_length_words for chord in self.chords
                if word[-1].is_composable_with(chord)]
            # filter on admissibility
            new_words = [word for word in new_words
                         if max(Counter([c.top_line_segment.knot_label for c in word]).values()) == 1]
            composable_admissible_words.append(new_words)
        cyclic_composable_admissible_words = [
            w for wl in range(len(composable_admissible_words))
            for w in composable_admissible_words[wl]
            if w[-1].is_composable_with(w[0])]
        self.rsft_generators = [WordOfChords(word) for word in cyclic_composable_admissible_words]
