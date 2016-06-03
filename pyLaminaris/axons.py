import numpy as np
import neuron

neuron.load_mechanisms('neuron')
from neuron import h


def _rotation(theta):
    R = np.zeros((3, 3))
    cx, cy, cz = np.cos(theta)
    sx, sy, sz = np.sin(theta)
    R.flat = (cx * cz - sx * cy * sz, cx * sz + sx * cy * cz, sx * sy,
              -sx * cz - cx * cy * sz, -sx * sz + cx * cy * cz,
              cx * sy, sy * sz, -sy * cz, cy)
    return R


def _gather_recursive(seg):
    if len(seg.children) == 0:
        return [seg]
    else:
        ret = [seg]
        for c in seg.children:
            ret.extend(_gather_recursive(c))
        return ret


class Segment:
    def __init__(self, parent=None, node_pitch=77.,
                 start=np.array([0., 0., 0.]), end=np.array([3 * 77., 0., 0.]), record=True):

        start = np.array(start)
        end = np.array(end)

        self.nnode = int(np.ceil(np.linalg.norm(end - start) / node_pitch))

        self.record = record

        if parent is None:
            self.order = 0
        else:
            self.order = parent.order + 1

        self.parent = parent
        self.children = []
        self.start = start
        self.end = end

        self.nodes = [h.Section() for x in range(self.nnode)]
        self.myelin = [h.Section() for x in range(self.nnode)]
        self.sections = []
        for n in range(self.nnode):
            self.sections.append(self.nodes[n])
            self.sections.append(self.myelin[n])
        # To connect two sections, call the connect() method of the child
        # Section object with the parent section as the argument:
        for n in range(self.nnode):
            self.nodes[n].connect(self.myelin[n])
        for n in range(self.nnode - 1):
            self.myelin[n + 1].connect(self.nodes[n])

        if self.parent is not None:
            self.myelin[0].connect(self.parent.nodes[-1])

        self._sections_setup()

        self.comps = np.array([comp for n in self.sections for comp in n.allseg()])
        self.comp_areas = np.array([c.area() for c in self.comps])
        self.comps = self.comps[self.comp_areas > 0]
        self.comp_areas = self.comp_areas[self.comp_areas > 0]
        self.ncomps = len(self.comps)

        if self.record:
            self.rec_i_mem = [h.Vector() for i in self.comps]
            [self.rec_i_mem[i].record(comp._ref_i_membrane)
             for i, comp in enumerate(self.comps)]

        self.node_locations = [self.start + (n + 1.) / self.ncomps * (self.end - self.start)
                               for n in range(self.ncomps)]

    def get_instantaneous_imem(self):
        return np.array([[self.nodes[i](0.5).i_membrane for i in range(self.nnode)]]).T

    def add_branch(self, direction=np.array([3 * 77., 0., 0.]),
                   angles=np.array([0., 0., 30.])):
        R1 = _rotation(angles / 360. * np.pi)  # half the angle in rad
        R2 = _rotation(-angles / 360. * np.pi)
        if direction.shape[0] == 3:
            direction = np.vstack((direction, direction))
        self.children.append(
            Segment(self, start=self.end, end=self.end + np.dot(direction[0, :], R1), record=self.record))
        self.children.append(
            Segment(self, start=self.end, end=self.end + np.dot(direction[1, :], R2), record=self.record))

    def _sections_setup(self):
        L_myelin = 75.
        L_node = 2.0
        diam_myelin = 1.3
        diam_node = 1.3

        area_node = np.pi * L_node * diam_node * 1.e-8  # in cm^2
        for sec in self.myelin:
            sec.nseg = 10
            sec.L = L_myelin
            sec.diam = diam_myelin

            sec.insert('leak')
            sec.insert('extracellular')

            sec.erev_leak = -20.
            sec.g_leak = 1e-6
            # sec.g_pas = 5.60e-9*l2a(sec.diam) # end up as S/cm2
            sec.cm = 1e-3  # 1.87e-11*l2a(sec.diam)*1e6 # end up as uF/cm2
            sec.Ra = 50.  # Ra(sec.diam)
        n = 0
        for sec in self.nodes:
            sec.nseg = 1
            sec.L = L_node
            sec.diam = diam_node

            sec.insert('leak')
            sec.insert('nahh')
            sec.insert('klva')
            sec.insert('khva')
            sec.insert('extracellular')
            if n == 0:
                sec.gnabar_nahh = 2.4
            else:
                sec.gnabar_nahh = .8
            # sec.gnabar_na = .8*(self.order+1)
            # sec.alphahVHalf_na = -58.
            # sec.betahVHalf_na = -58.
            # sec.alphahK_na = 5.
            #sec.betah_na = 5.
            #sec.alphamVHalf_na = -39.
            #sec.betamVHalf_na = -39.
            #sec.alphamK_na = 5.
            #sec.betamK_na = 5.

            sec.gkbar_klva = 0.1
            sec.alphaVHalf_klva = -50.
            sec.betaVHalf_klva = -50.
            sec.q10_klva = 3.0

            sec.gkbar_khva = 1.5
            sec.q10_khva = 3.0

            sec.g_leak = 1e-3
            sec.erev_leak = -20.

            sec.ena = 50.  #115 - 65
            sec.ek = -80.  #-12 - 65
            sec.cm = 1.0  #3.14e-9*l2a(sec.diam)*1e6 # end up as uF/cm2
            sec.Ra = 50.  #Ra(sec.diam)


class Tree:
    def __init__(self, depth=5, root_point=np.array([3000., 0., 0.]), record=True, mtype='long'):
        self.mtype = mtype
        self.record = record
        #self.root = Segment(end=root_point, record=record)
        self.root = Segment(start=[np.random.random()*77., root_point[1], root_point[2]],
                            end=root_point, record=self.record)
        self.build_deterministic(depth)
        self.virgin = True
        self.delay = [1.]

    def build_deterministic(self, depth):
        def gamma_m_v(m,st):
            v = float(st)**2
            t = v/float(m)
            k = float(m)/t
            return np.random.gamma(k,t)
        while True:
            to_extend = [s for s in self.leaves() if s.order < depth]
            if len(to_extend) == 0:
                break
            for s in to_extend:
                if self.mtype == 'pop_sync' and depth == 1:
                    s.add_branch(direction=np.array(
                        [[np.random.gamma(20,25), 0., 0.],
                         [np.random.gamma(20,25), 0., 0.]]),
                        angles=np.array([0., 0., 0.]))
                elif self.mtype == 'pop_sync' and depth > 1:
                    if s.order == depth-1:
                        s.add_branch( direction=np.array([
                            [gamma_m_v(400,10), 0., 0.],
                            [gamma_m_v(400,10), 0., 0.]]
                            ),
                                angles=np.array([0., 0., 0.]))
                    else:
                        s.add_branch(direction=np.array(
                            [[gamma_m_v(400,300), 0., 0.],
                             [gamma_m_v(400,300), 0., 0.]]),
                            angles=np.array([0., 0., 0.]))
                elif self.mtype == 'bif':
                    if s.order == depth-1:
                        s.add_branch( direction=np.array([
                            [10000., 0., 0.],
                            [10000., 0., 0.]]
                            ),
                                angles=np.array([0., 0., 0.]))
                    else:
                        s.add_branch(direction=np.array(
                            [[100., 0., 0.],
                             [100., 0., 0.]]),
                            angles=np.array([0., 0., 0.]))
                elif self.mtype == 'bif_term':
                    if s.order == depth-1:
                        s.add_branch( direction=np.array([
                            [700., 0., 0.],
                            [700., 0., 0.]]
                            ),
                                angles=np.array([0., 0., 0.]))
                    else:
                        s.add_branch(direction=np.array(
                            [[100., 0., 0.],
                             [100., 0., 0.]]),
                            angles=np.array([0., 0., 0.]))
                else:
                    s.add_branch( direction=np.array([10000., 0., 0.]),
                            angles=np.array([0., 0., 0.]))

    def segments(self):
        return _gather_recursive(self.root)

    def spatial_layout(self):
        bif = []
        cont = []
        term = []

        for s in self.segments():
            if len(s.children) == 0:
                term.append(s.end)
            else:
                bif.append(s.end)

            cont.extend(s.node_locations)

        return np.array(bif), np.array(cont), np.array(term)

    def leaves(self):
        return [s for s in self.segments() if len(s.children) == 0]

    def nodes_imem_loc(self, with_root=True):
        import numpy as np

        segs = self.segments()
        imem = []
        loc = []
        # exclude first sections of root segment
        #if segs[0].record:
        #    imem.extend(segs[0].rec_i_mem[10:])
        #else:
        #    imem.extend(segs[0].get_instantaneous_imem()[10:])
        #loc.extend(segs[0].node_locations[10:])
        if with_root:
            start = 0
        else:
            start = 1
        for s in segs[start:]:
            # TODO: here is where we need to differentiate record or not.
            # FIXME: Should be failing a test!!!!
            if s.record:
                imem.extend(s.rec_i_mem * s.comp_areas[:, np.newaxis])
            else:
                imem.extend(s.get_instantaneous_imem())
            loc.extend(s.node_locations)
        imem = [np.array(i) for i in imem]

        return np.array(imem), np.array(loc)

    def draw_2d(self, alpha=1.0, color='b'):
        from matplotlib import pyplot as plt
        import numpy as np

        for s in self.segments():
            l = np.vstack((s.start, s.end))
            plt.plot(l[:, 0], l[:, 1], 'k', alpha=alpha, color=color)

            n = np.vstack(s.node_locations)
            plt.plot(n[:, 0], n[:, 1], 'k.', alpha=alpha, color=color)

        plt.gca().set_aspect('equal')


class ProbTree(Tree):
    def __init__(self, **kwargs):
        from scipy import stats

        if 'record' in kwargs.keys():
            self.record = kwargs['record']
            kwargs.pop('record')
        else:
            self.record = True

        if 'structure' in kwargs.keys():
            structure = kwargs['structure']
            kwargs.pop('structure')
        else:
            structure = 'prob'

        if 'mode' in kwargs.keys():
            mode = kwargs['mode']
            kwargs.pop('mode')
        else:
            mode = 'pulse'

        if 'root_point' in kwargs.keys():
            if type(kwargs['root_point']) is float or type(kwargs['root_point']) is int:
                self.root = Segment(end=np.array([kwargs['root_point'], 0, 0]), record=self.record)
            else:
                self.root = Segment(start=[np.random.random()*77., kwargs['root_point'][1], kwargs['root_point'][2]],
                                    end=kwargs['root_point'], record=self.record)
            kwargs.pop('root_point')
        else:
            self.root = Segment(record=self.record)

        if structure == 'prob':
            self.build_prob(**kwargs)
        elif structure == 'logistic':
            self.build_logistic(**kwargs['logistic_params'])

        self.virgin = True

    def build_logistic(self,
                       bif_center=15000., bif_sigma=200., bif_amp=500.,
                       ter_center=15500., ter_sigma=100.,
                       ang_mean=20., ang_var=5.
                       ):

        def bif(x, A, mu, sigma):
            return A / sigma * np.exp((x - mu) / sigma) / (1 + np.exp((x - mu) / sigma)) ** 2

        def term(x, mu, sigma):
            return 1. - 1. / (1 + np.exp((x - mu) / sigma))

        def build_recursive(node):
            root_x = node.end[0]
            p_bif = bif(root_x, bif_amp, bif_center, bif_sigma)
            p_term = term(root_x, ter_center, ter_sigma)
            if np.random.rand() < p_bif / (p_bif + p_term):
                dir_1 = 77. * 6 * np.random.rand()
                dir_2 = 77. * 6 * np.random.rand()
                node.add_branch(
                    direction=np.array([[dir_1, 0., 0.], [dir_2, 0., 0.]]),
                    angles=np.array([0., 0., ang_mean + np.random.rand() * ang_var]))
                for s in node.children:
                    build_recursive(s)

        build_recursive(self.root)

    def build_prob(self, depth=5, p=[1., .9, .7, .5, .5, .4],
                   dir_mean=154., dir_var=77.,
                   ang_mean=25, ang_var=20):
        # ang_mean=37.5, ang_var=12.5):
        def build_recursive(node, depth, p, dir_mean, dir_var, ang_mean, ang_var):
            if (node.order < depth) and (np.random.rand() < p[node.order]):
                node.add_branch(
                    direction=np.array([dir_mean + np.random.rand() * dir_var, 0., 0.]),
                    angles=np.array([0., 0., ang_mean + np.random.rand() * ang_var]))
                for s in node.children:
                    build_recursive(s, depth, p, dir_mean, dir_var, ang_mean, ang_var)

        build_recursive(self.root, depth, p, dir_mean, dir_var, ang_mean, ang_var)
