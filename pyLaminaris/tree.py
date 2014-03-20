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


def _gather_rec(seg):
    if len(seg.children) == 0:
        return [seg]
    else:
        ret = [seg]
        for c in seg.children:
            ret.extend(_gather_rec(c))
        return ret


class Segment:
    def __init__(self, parent=None, node_pitch=77.,
                 start=np.array([0., 0., 0.]), end=np.array([3 * 77., 0., 0.])):

        nnode = int(np.ceil(np.linalg.norm(end - start) / node_pitch))

        if parent is None:
            self.order = 0
        else:
            self.order = parent.order + 1

        self.parent = parent
        self.children = []
        self.start = start
        self.end = end

        self.nodes = [h.Section() for x in range(nnode)]
        self.myelin = [h.Section() for x in range(nnode)]
        # To connect two sections, call the connect() method of the child
        # Section object with the parent section as the argument:
        for n in range(nnode):
            self.nodes[n].connect(self.myelin[n])
        for n in range(nnode - 1):
            self.myelin[n + 1].connect(self.nodes[n])

        if self.parent is not None:
            self.myelin[0].connect(self.parent.nodes[-1])

        self._sections_setup()

        self.rec_i_mem = [h.Vector() for i in range(nnode)]
        [self.rec_i_mem[i].record(self.nodes[i](0.5)._ref_i_membrane)
         for i in range(nnode)]

        self.node_locations = [self.start + (n + 1.) / nnode * (self.end - self.start)
                               for n in range(nnode)]


    def add_branch(self, direction=np.array([3 * 77., 0., 0.]),
                   angles=np.array([0., 0., 30.])):
        R1 = _rotation(angles / 360. * np.pi)  #half the angle in rad
        R2 = _rotation(-angles / 360. * np.pi)
        if direction.shape[0] == 3:
            direction = np.vstack((direction, direction))
        self.children.append(
            Segment(self, start=self.end, end=self.end + np.dot(direction[0, :], R1)))
        self.children.append(
            Segment(self, start=self.end, end=self.end + np.dot(direction[1, :], R2)))

    def _sections_setup(self):
        L_myelin = 75.
        L_node = 2.0
        diam_myelin = 1.3
        diam_node = 1.3

        area_node = np.pi * L_node * diam_node * 1.e-8  #in cm^2
        for sec in self.myelin:
            sec.nseg = 10
            sec.L = L_myelin
            sec.diam = diam_myelin

            sec.insert('leak')
            sec.insert('extracellular')

            sec.erev_leak = -20.
            sec.g_leak = 1.0e-6
            #sec.g_pas = 5.60e-9*l2a(sec.diam) # end up as S/cm2
            sec.cm = 0.001  #1.87e-11*l2a(sec.diam)*1e6 # end up as uF/cm2
            sec.Ra = 50.  #Ra(sec.diam)

        for sec in self.nodes:
            sec.nseg = 1
            sec.L = L_node
            sec.diam = diam_node

            sec.insert('leak')
            sec.insert('na')
            sec.insert('klva')
            sec.insert('khva')
            sec.insert('extracellular')

            sec.gnabar_na = .8
            #sec.alphahVHalf_na = -58.
            #sec.betahVHalf_na = -58.
            #sec.alphahK_na = 5.
            #sec.betahK_na = 5.
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

            sec.g_leak = 0.0006
            sec.erev_leak = -20.

            sec.ena = 50.  #115 - 65
            sec.ek = -80.  #-12 - 65
            sec.cm = 1.0  #3.14e-9*l2a(sec.diam)*1e6 # end up as uF/cm2
            sec.Ra = 50.  #Ra(sec.diam)


class Tree:
    def __init__(self, depth=5, root_point=np.array([3000., 0., 0.])):
        self.root = Segment(end=root_point)
        self.build_deterministic(depth)
        self.virgin = True
        self.delay = [1.]

    def build_deterministic(self, depth):
        while True:
            to_extend = [s for s in self.leaves() if s.order < depth]
            if len(to_extend) == 0:
                break
            for s in to_extend:
                s.add_branch()

    def segments(self):
        return _gather_rec(self.root)

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

    def nodes_imem_loc(self):
        import numpy as np

        segs = self.segments()
        imem = []
        loc = []
        for s in segs:
            imem.extend(s.rec_i_mem)
            loc.extend(s.node_locations)
        imem = [np.array(i) for i in imem]

        return np.array(imem), np.array(loc)

    def run_sim(self):
        if not self.virgin:
            print "warning, rerunning simulation on non-virgin tree"
        stims = []
        for d in self.delay:
            stim = h.IClamp(self.root.nodes[0](0.5))
            # Setting up stimulation
            stim.delay = d
            stim.amp = 1.
            stim.dur = 0.1
            stims.append(stim)
        h.celsius = 40.
        h.dt = 0.0025
        h.finitialize(-75)
        neuron.init()
        neuron.run(20)
        self.virgin = False

    def draw_2d(self):
        from matplotlib import pyplot as plt
        import numpy as np

        for s in self.segments():
            l = np.vstack((s.start, s.end))
            plt.plot(l[:, 0], l[:, 1], 'k')

            n = np.vstack(s.node_locations)
            plt.plot(n[:, 0], n[:, 1], 'k.')

        plt.gca().set_aspect('equal')

    def calc_fields(self, x=3000., y=400., z=0.):
        assert not self.virgin, "need to run simulation first"
        import calc_lfp
        import numpy as np

        imem, loc = self.nodes_imem_loc()
        if type(x) is float:
            x = np.array([x])
        sample_loc = np.vstack([x,
                                y * np.ones(x.shape),
                                z * np.ones(x.shape)]).T

        return calc_lfp.calc_kernels_3d(loc, imem, sample_loc)


class ProbTree(Tree):
    def __init__(self, **kwargs):
        from scipy import stats

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
            self.root = Segment(end=np.array([kwargs['root_point'], 0, 0]))
            kwargs.pop('root_point')
        else:
            self.root = Segment()

        if structure == 'prob':
            self.build_prob(**kwargs)
        elif structure == 'logistic':
            self.build_logistic(**kwargs)

        self.virgin = True
        #self.delay = 1.5 + 0.25*np.random.randn()

        if mode == 'pulse':
            class rv(stats.rv_continuous):
                def _pdf(self, x):
                    return (1. + np.sin(6.28 * 3.33 * x)) * np.exp(-(x * x / (2 * .5 * .5))) / 1.253

            myrv = rv()
            self.delay = [0.]
            while self.delay[0] < 0.2:
                self.delay[0] = myrv.rvs(1.5)
        elif mode == 'spont':
            import helper

            self.delay = helper.hom_poisson(0.1, 100, .1)
        elif mode == 'click':
            import helper

            def click(t):
                return 0.1 + 0.9 * np.exp(np.sin(6 * np.pi * t) - 1. - (t - 35) ** 2)

            self.delay = helper.inhom_poisson(click, 50., 0., 200.)
        elif mode == 'mod_click':
            import helper

            def mod_click(t):
                return (0.1
                        + 2.9 * np.exp(np.sin(4 * 2 * np.pi * t) - 1. - 16 * (t - 10) ** 2)
                        + 1.5 * np.exp(np.sin(4 * 2 * np.pi * t) - 1. - 16 * (t - 11) ** 2)
                        - 0.1 * np.exp(-16 * (t - 10.5) ** 2))

            self.delay = helper.inhom_poisson(mod_click, 20., 0., 6.)

    def build_logistic(self,
                       bif_center=5000., bif_sigma=200., bif_amp=500.,
                       ter_center=5500., ter_sigma=100.,
                       ang_mean=20., ang_var=0.
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
                   ang_mean=37.5, ang_var=12.5):
        def build_recursive(node, depth, p, dir_mean, dir_var, ang_mean, ang_var):
            if (node.order < depth) and (np.random.rand() < p[node.order]):
                node.add_branch(
                    direction=np.array([dir_mean + np.random.rand() * dir_var, 0., 0.]),
                    angles=np.array([0., 0., ang_mean + np.random.rand() * ang_var]))
                for s in node.children:
                    build_recursive(s, depth, p, dir_mean, dir_var, ang_mean, ang_var)

        build_recursive(self.root, depth, p, dir_mean, dir_var, ang_mean, ang_var)
