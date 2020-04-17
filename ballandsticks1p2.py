# =============================================================================
# BALL AND STICKS CLASSES 1 (1 RGC INPUT AND TRIAD) - PARAMETERISATION 2
# -----------------------------------------------------------------------------
# Patrick McCarthy, Laura Lazzari and Jonathan Martin, March 2020 
# ----------------------------------------------------------------------------- 
# This file defines classes for an LGN interneuron with 1 dendrite and a 
# 1 axon, to be used in conjunction with the file 'model1.py' to simulate a
# small mouse LGN network with 1 RCG input to a relay cell, which is
# inhibited by the interneuron via axosomatic (F1) and triadic (F2) inhibitory
# synapses. This is parameterisation 2 of the model classes.
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
        self.dend_p = h.Section(name='proximal_part_of_dendrite',cell=self)
        self.dend_d = h.Section(name='distal_part_of_dendrite',cell=self)
        self.all = [self.soma, self.axon_p, self.axon_d, self.dend_p, self.dend_d]
        self.exc_soma = [self.axon_p, self.axon_d, self.dend_p, self.dend_d]
        self.axon_p.connect(self.soma(0))
        self.axon_d.connect(self.axon_p(1))
        self.dend_p.connect(self.soma(1)) 
        self.dend_d.connect(self.dend_p(1))
        self.soma.nseg = 11
        self.soma.L = 15.3 
        self.soma.diam = 17.4 
        self.axon_p.nseg = 11
        self.axon_p.L = 100
        taper_diam(self.axon_p,4,0.3)
        self.axon_d.nseg = 11
        self.axon_d.L = 400 
        self.axon_d.diam = 0.3 
        self.dend_p.nseg = 11
        self.dend_p.L = 100
        taper_diam(self.dend_p,4,0.3)
        self.dend_d.nseg = 11 
        self.dend_d.L = 400 
        self.dend_d.diam = 0.3 
    
    # biophysics
    def _setup_biophysics(self):
        for sec in self.all: 
            sec.Ra = 113
            sec.cm = 1.1
        self.soma.insert('hh') # Hodgkin-Huxley kinetics
        for seg in self.soma:                       
            seg.hh.gnabar = 0.5  
            seg.hh.gkbar = 0.1  
            seg.hh.gl = 0.0003       
            seg.hh.el = (-54.3*mV) # Nernst potential
        self.soma.insert('pas')
        for seg in self.soma:
            seg.pas.g = 0.0001 
            seg.pas.e = -60*mV    
         
        for sec in [self.dend_p,self.dend_d,self.axon_p,self.axon_d]:
            sec.insert('pas') # Passive membrane properties
            sec.insert('hh')
            for seg in sec:
                seg.pas.g = 0.005  
                seg.pas.e = (-65*mV) 
                seg.hh.gnabar = 0.5  # arbitrarily chosen
                seg.hh.gkbar = 0.1  
                seg.hh.gl = 0.0003       
                seg.hh.el = (-54.3*mV) # 
                
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
            seg.hh.el = -54.3*mV
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
            sec.ena=50
            sec.ek=-90
            sec.insert('ican')  # CAN-channel from Zhu et al. 99a

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