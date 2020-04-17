# =============================================================================
# BALL AND STICKS CLASSES 2 (3 RGC INPUTS AND TRIADS)
# -----------------------------------------------------------------------------
# Patrick McCarthy, Laura Lazzari and Jonathan Martin, March 2020 
# ----------------------------------------------------------------------------- 
# This file defines classes for an LGN interneuron with 3 dendrites and a 
# 1 axon, to be used in conjunction with the file 'model2.py' to simulate a
# small mouse LGN network with 3 RCG inputs to a relay cell, which is
# inhibited by the interneuron via axosomatic (F1) and triadic (F2) inhibitory
# synapses.
# -----------------------------------------------------------------------------
# We thank Prof. Gaute Einevoll, Dr. Thomas Heiberg and Dr. Geir Halnes for
# providing us with the code on which we based this model.
# =============================================================================

from neuron import h
from neuron.units import mV
import LFPy
h.load_file('stdrun.hoc')

class Interneuron(LFPy.TemplateCell):
    
    # constructor 
    def __init__(self, gid, x, y, z, theta):
        self._gid = gid 
        self._setup_morphology()
        self._setup_biophysics()
        self.x = self.y = self.z = 0
        h.define_shape()
        self._rotate_z(theta)                             
        self._set_position(x, y, z)                     
        
    # morphology
    def _setup_morphology(self): 
        self.soma = h.Section(name='soma', cell=self)
        self.axon_p = h.Section(name='proximal_part_of_axon',cell=self)
        self.axon_d = h.Section(name='distal_part_of_axon',cell=self)
        self.dend1_p = h.Section(name='proximal_part_of_dendrite 1',cell=self)
        self.dend1_d = h.Section(name='distal_part_of_dendrite 1',cell=self)
        self.dend2_p = h.Section(name='proximal_part_of_dendrite 2',cell=self)
        self.dend2_d = h.Section(name='distal_part_of_dendrite 2',cell=self)
        self.dend3_p = h.Section(name='proximal_part_of_dendrite 3',cell=self)
        self.dend3_d = h.Section(name='distal_part_of_dendrite 3',cell=self)
        self.all = [self.soma, self.axon_p, self.axon_d, self.dend1_p, self.dend1_d, self.dend2_p, self.dend2_d, self.dend3_p, self.dend3_d]
        self.exc_soma = [self.axon_p, self.axon_p, self.axon_d, self.dend1_p, self.dend1_d, self.dend2_p, self.dend2_d, self.dend3_p, self.dend3_d]
        self.axon_p.connect(self.soma(0))
        self.axon_d.connect(self.axon_p(1))
        self.dend1_p.connect(self.soma(0.3)) 
        self.dend1_d.connect(self.dend1_p(1))
        self.dend2_p.connect(self.soma(0.6)) 
        self.dend2_d.connect(self.dend2_p(1))
        self.dend3_p.connect(self.soma(0.9)) 
        self.dend3_d.connect(self.dend2_p(1))
        self.soma.nseg = 11
        self.soma.L = 15.3 
        self.soma.diam = 17.4 
        self.axon_p.nseg = 11
        self.axon_p.L = 100
        taper_diam(self.axon_p,4,0.3)
        self.axon_d.nseg = 11
        self.axon_d.L = 400 
        self.axon_d.diam = 0.3 
        self.dend1_p.nseg = 11
        self.dend1_p.L = 100
        taper_diam(self.dend1_p,4,0.3)
        self.dend1_d.nseg = 11 
        self.dend1_d.L = 400 
        self.dend1_d.diam = 0.3 
        self.dend2_p.nseg = 11
        self.dend2_p.L = 100
        taper_diam(self.dend2_p,4,0.3)
        self.dend2_d.nseg = 11 
        self.dend2_d.L = 400 
        self.dend2_d.diam = 0.3 
        self.dend3_p.nseg = 11
        self.dend3_p.L = 100
        taper_diam(self.dend3_p,4,0.3)
        self.dend3_d.nseg = 11 
        self.dend3_d.L = 400 
        self.dend3_d.diam = 0.3 
    
    # biophysics
    def _setup_biophysics(self):
        for sec in self.all: 
            sec.Ra = 250
            sec.cm = 1
        self.soma.insert('hh') # Hodgkin-Huxley kinetics
        for seg in self.soma:                       
            seg.hh.gnabar = 0.05  
            seg.hh.gkbar = 0.05  
            seg.hh.gl = 0.0003       
            seg.hh.el = (-50*mV) # Nernst potential
        self.soma.insert('pas')
        for seg in self.soma:
            seg.pas.g = 0.0001 
            seg.pas.e = -60*mV    
         
        for sec in self.exc_soma:
            sec.insert('pas') # Passive membrane properties
            sec.insert('hh')
            for seg in sec:
                seg.pas.g = 0.005  
                seg.pas.e = (-65*mV) 
                seg.hh.gnabar = 0.5  # arbitrarily chosen
                seg.hh.gkbar = 0.1  
                seg.hh.gl = 0.0003       
                seg.hh.el = (-50*mV) # 
                
        for sec in self.all: # add rest of ion channel types
            sec.insert('iar')
            sec.insert('Cad')   # Calsium pool, Zhu et al.
            sec.insert('ical')  # L-type Ca-current, using pool in Cad
            sec.insert('it2')   # t-type Ca- current, using pool in Cad
            sec.insert('iahp')  # potassium current, slow, Ca-dependent, Zhu et al.
            sec.insert('hh2')
            sec.ena=50 # Sodium reversal potential
            sec.ek=-90 # Potassium reversal potential
            sec.insert('ican')  # CAN-channel from Zhu et al. 99a

        
    # specify how sections are to be diplayed
    def __repr__(self):
        return 'Interneuron[{}]'.format(self._gid)
    
    # position in 3D space
    def _set_position(self, x, y, z):
        for sec in self.all: 
            for i in range(sec.n3d()): 
                sec.pt3dchange(i,
                               x - self.x + sec.x3d(i),
                               y - self.y + sec.y3d(i),
                               z - self.z + sec.z3d(i),
                              sec.diam3d(i))
        self.x, self.y, self.z = x, y, z
        
    # rotate cell about z-axis
    def _rotate_z(self, theta):
        for sec in self.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = h.cos(theta)
                s = h.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))
                
class RelayCell(LFPy.TemplateCell):
    
    # constructor 
    def __init__(self, gid, x, y, z, theta):
        self._gid = gid
        self._setup_morphology()
        self._setup_biophysics()
        self.x = self.y = self.z = 0
        h.define_shape()
        self._rotate_z(theta)                             
        self._set_position(x, y, z)                     
        
        
    # morphology
    def _setup_morphology(self): 
        self.soma = h.Section(name='soma', cell=self)
        self.all = [self.soma] 
        self.soma.L = 50 
        self.soma.diam = 47 
        self.soma.nseg = 11 
    
    # biophysics
    def _setup_biophysics(self):
        for sec in self.all: 
            sec.Ra = 1000
            sec.cm = 1 
        self.soma.insert('hh')                                   
        for seg in self.soma:                            
            seg.hh.gnabar = 0.1
            seg.hh.gkbar = 0.025  
            seg.hh.gl = 0.0001
            seg.hh.el = -50*mV
        self.soma.insert('pas')
        for seg in self.soma:
            seg.pas.g = 0.0001 
            seg.pas.e = -60*mV    
            
        for sec in self.all: # add rest of ion channel types
            sec.insert('iar')
            sec.insert('Cad')   # Calsium pool, Zhu et al.
            sec.insert('ical')  # L-type Ca-current, using pool in Cad
            sec.insert('it2')   # t-type Ca- current, using pool in Cad
            sec.insert('iahp')  # potassium current, slow, Ca-dependent, Zhu et al.
            sec.insert('hh2')
            sec.insert('ican')  # CAN-channel from Zhu et al. 99a
            sec.ena=50
            sec.ek=-90

    # specify how sections are to be diplayed
    def __repr__(self):
        return 'RelayCell[{}]'.format(self._gid)
    
    # specify position in 3D space
    def _set_position(self, x, y, z):
        for sec in self.all: 
            for i in range(sec.n3d()): 
                sec.pt3dchange(i,
                               x - self.x + sec.x3d(i),
                               y - self.y + sec.y3d(i),
                               z - self.z + sec.z3d(i),
                              sec.diam3d(i))
        self.x, self.y, self.z = x, y, z
        
    # rotate cell about z-axis
    def _rotate_z(self, theta):
        for sec in self.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = h.cos(theta)
                s = h.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))

     
def taper_diam(sec,zero_bound,one_bound):
    dx=1.0/(sec.nseg)
    x=dx/2
    for seg in sec:
        seg.diam=(one_bound-zero_bound)*x+zero_bound
        x+=dx