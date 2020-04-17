# =============================================================================
# MODEL 1 (1 RGC INPUT AND TRIAD)
# -----------------------------------------------------------------------------
# Patrick McCarthy, Laura Lazzari and Jonathan Martin, March 2020
# ----------------------------------------------------------------------------- 
# This file is to be used in conjunction with either 'ballandsticks1p1.py'
# or 'ballandsticks1p2.py'. It implements classes from these files to create
# cell objects and then adds synapses to simulate both the input stimuli and 
# excitatory and inhibitory connections between cells.
# -----------------------------------------------------------------------------
# We thank Prof. Gaute Einevoll, Dr. Thomas Heiberg and Dr. Geir Halnes for
# providing us with the code on which we based this model.
# =============================================================================

# import libraries
import ballandsticks1p1 as bs # make sure you have the 'ballandsticks.py' file in the same directory or else this won't work and the code won't run
from neuron import h
from neuron.units import ms, mV

h.load_file('stdrun.hoc')

# objects for interneuron and relay cell
interneuron = bs.Interneuron(0,0,0,1,0) 
relaycell = bs.RelayCell(0,0,0,1,0) 

# vectors to store synapse and connections
syns = []
netcons = []

# stimulator object
stim = h.NetStim() 
stim.number = 1
stim.start = 5 * ms 
stim.interval = 4.5 * ms

# excitatory synapse on relay cell (part 1 of triad) 1
rc_exc = h.Exp2Syn(relaycell.soma(0.95))
rc_exc_con = h.NetCon(stim, rc_exc) 
rc_exc_con.delay = 0 * ms
rc_exc_con.weight[0] = 5 # Hieberg: 11.6
rc_exc.e = 42 * mV
rc_exc.tau1 = 1 * ms 
rc_exc.tau2 = 2 * ms 
syns.append(rc_exc)
netcons.append(rc_exc_con)

# proximal excitatory axodendritic synapse on interneuron
in_exc = h.Exp2Syn(interneuron.dend_p(0.1)) 
in_exc_con = h.NetCon(stim, in_exc) 
in_exc_con.weight[0] = 0.6
in_exc_con.delay = 0 * ms 
in_exc.e = 42 * mV
in_exc.tau1 = 1.6 * ms 
in_exc.tau2 = 3.6 * ms 
syns.append(in_exc)
netcons.append(in_exc_con)

# distal excitatory axodendritic synapse on interneuron (part 2 of triad)
triad1 = h.Exp2Syn(interneuron.dend_d(1)) 
triad1_con = h.NetCon(stim, triad1) 
triad1_con.weight[0] = 2
triad1_con.delay = 0 * ms 
triad1.e = 42 * mV
triad1.tau1 = 0.3 * ms 
triad1.tau2 = 2 * ms 
syns.append(triad1)
netcons.append(triad1_con)

# inhibitory dendrodendritic synapse on relay cell (part 3 of triad)
triad2 = h.Exp2Syn(relaycell.soma(0.9)) 
triad2_con = h.NetCon(interneuron.dend_d(0.9)._ref_v, triad2, sec=interneuron.dend_d) 
triad2_con.weight[0] = 10
triad2_con.delay = 0.5 * ms 
triad2.e = -75 * mV
triad2.i = 10000
triad2.tau2 = 0.7 * ms 
triad2.tau2 = 4.2 * ms 
syns.append(triad2)
netcons.append(triad2_con)

# inhibitory axosomatic synapse between interneuron and relay cell
rc_inh = h.Exp2Syn(relaycell.soma(0.1)) 
rc_inh_con = h.NetCon(interneuron.axon_d(1)._ref_v, rc_inh, sec=interneuron.axon_d) 
rc_inh_con.weight[0] = 10
rc_inh_con.delay = 1  * ms 
rc_inh.e = -75 * mV
rc_inh.i = 10000
rc_inh.tau1 = 0.7 * ms 
rc_inh.tau2 = 4.2 * ms 
syns.append(rc_inh)
netcons.append(rc_inh_con)
