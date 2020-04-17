# =============================================================================
# MODEL 4(3 RGC INPUTS AND TRIADS)
# -----------------------------------------------------------------------------
# Patrick McCarthy, Laura Lazzari and Jonathan Martin, March 2020 
# ----------------------------------------------------------------------------- 
# This file is to be used in conjunction with 'ballandsticks2.py'. It 
# implements classes from this file to create cell objects and then adds 
# AlphaSynapse objects to simulate input stimuli and excitatory with NetCon objects 
# as inhibitory connections between cells.
# -----------------------------------------------------------------------------
# We thank Prof. Gaute Einevoll, Dr. Thomas Heiberg and Dr. Geir Halnes for
# providing us with the code on which we based this model.
# =============================================================================

# IMPORT LIBRARIES
import ballandsticks2 as bs2 # make sure you have the 'ballandsticks.py' file in the same directory or else this won't work and the code won't run
from neuron import h
from neuron.units import ms, mV
import matplotlib.pyplot as plt


h.load_file('stdrun.hoc')

# objects for interneuron and relay cell
interneuron = bs2.Interneuron(0,0,0,1,0) 
relaycell = bs2.RelayCell(0,0,0,1,0) 

# vectors to store synapse and connections
syns = []
netcons = []

# STIMULATOR OBJECTS

## stimulator 1
#stim1 = h.NetStim() 
#stim1.number = 1
#stim1.start = 5 * ms 
#stim1.interval = 0.5 * ms
#
## stimulator 2
#stim2 = h.NetStim() 
#stim2.number = 1
#stim2.start = 5 * ms 
#stim2.interval = 0.5 * ms
#
## stimulator 2
#stim3 = h.NetStim() 
#stim3.number = 1
#stim3.start = 5 * ms 
#stim3.interval = 0.5 * ms

# SYNAPSES

# set 1

# excitatory synapse on relay cell - technically part of triad but named differently for clarity
rc_exc1 = h.AlphaSynapse(relaycell.soma(0.28))
rc_exc1.onset = 5 * ms
rc_exc1.tau = 0.2 * ms
rc_exc1.gmax = 5
rc_exc1.i = 12
syns.append(rc_exc1)

# proximal excitatory axodendritic synapse on interneuron
in_exc1 = h.AlphaSynapse(interneuron.dend1_p(0.1)) 
in_exc1.onset = 5 * ms
in_exc1.tau = 0.2 * ms
in_exc1.gmax = 5
in_exc1.i = 12
syns.append(in_exc1)

# distal excitatory axodendritic synapse on interneuron 
triad1_1 = h.AlphaSynapse(interneuron.dend1_d(1)) 
triad1_1.onset = 5 * ms
triad1_1.tau = 0.5 * ms
triad1_1.gmax = 5
triad1_1.i = 12
syns.append(triad1_1)

# inhibitory dendrodendritic synapse on relay cell
triad1_2 = h.Exp2Syn(relaycell.soma(0.3)) 
triad1_2_con = h.NetCon(interneuron.dend1_d(0.99)._ref_v, triad1_2, sec=interneuron.dend1_d) 
triad1_2_con.weight[0] = 1.5
triad1_2_con.delay = 0.5 * ms 
triad1_2.e = -75 * mV
triad1_2.tau2 = 0.7 * ms 
triad1_2.tau2 = 4.2 * ms 
syns.append(triad1_2)
netcons.append(triad1_2_con)

# set 2

# excitatory synapse on relay cell - technically part of triad but named differently for clarity
rc_exc2 = h.AlphaSynapse(relaycell.soma(0.58))
rc_exc2.onset = 5 * ms
rc_exc2.tau = 0.2 * ms
rc_exc2.gmax = 10
rc_exc2.i = 12
syns.append(rc_exc2)

# proximal excitatory axodendritic synapse on interneuron
in_exc2 = h.AlphaSynapse(interneuron.dend2_p(0.1)) 
in_exc2.onset = 5 * ms
in_exc2.tau = 0.2 * ms
in_exc2.gmax = 10
in_exc2.i = 12
syns.append(in_exc2)

# distal excitatory axodendritic synapse on interneuron 
triad2_1 = h.AlphaSynapse(interneuron.dend2_d(1)) 
triad2_1.onset = 5 * ms
triad2_1.tau = 0.2 * ms
triad2_1.gmax = 10
triad2_1.i = 12
syns.append(triad2_1)

# inhibitory dendrodendritic synapse on relay cell 
triad2_2 = h.Exp2Syn(relaycell.soma(0.6)) 
triad2_2_con = h.NetCon(interneuron.dend2_d(0.99)._ref_v, triad2_2, sec=interneuron.dend2_d) 
triad2_2_con.weight[0] = 1.5
triad2_2_con.delay = 0.5 * ms 
triad2_2.e = -75 * mV
triad2_2.tau2 = 0.7 * ms 
triad2_2.tau2 = 4.2 * ms 
syns.append(triad2_2)
netcons.append(triad2_2_con)

# set 3

# excitatory synapse on relay cell - technically part of triad but named differently for clarity
rc_exc3 = h.AlphaSynapse(relaycell.soma(0.98))
rc_exc3.onset = 5 * ms
rc_exc3.tau = 0.2 * ms
rc_exc3.gmax = 15
rc_exc3.i = 12
syns.append(rc_exc3)

# proximal excitatory axodendritic synapse on interneuron
in_exc3 = h.AlphaSynapse(interneuron.dend3_p(0.1)) 
in_exc3.onset = 5 * ms
in_exc3.tau = 0.2 * ms
in_exc3.gmax = 15
in_exc3.i = 12
syns.append(in_exc3)

# distal excitatory axodendritic synapse on interneuron 
triad3_1 = h.AlphaSynapse(interneuron.dend3_d(1)) 
triad3_1.onset = 5 * ms
triad3_1.tau = 0.2 * ms
triad3_1.gmax = 15
triad3_1.i = 12
syns.append(triad3_1)

# inhibitory dendrodendritic synapse on relay cell 
triad3_2 = h.Exp2Syn(relaycell.soma(0.9)) 
triad3_2_con = h.NetCon(interneuron.dend3_d(0.99)._ref_v, triad3_2, sec=interneuron.dend3_d) 
triad3_2_con.weight[0] = 1.5
triad3_2_con.delay = 0.5 * ms 
triad3_2.e = -75 * mV
triad3_2.tau2 = 0.7 * ms 
triad3_2.tau2 = 4.2 * ms 
syns.append(triad3_2)
netcons.append(triad3_2_con)

# inhibitory axosomatic synapse between interneuron and relay cell
rc_inh = h.Exp2Syn(relaycell.soma(0.1)) 
rc_inh_con = h.NetCon(interneuron.axon_d(1)._ref_v, rc_inh, sec=interneuron.axon_d) 
rc_inh_con.weight[0] = 1.5
rc_inh_con.delay = 1  * ms 
rc_inh.e = -75 * mV
rc_inh.tau1 = 0.7 * ms 
rc_inh.tau2 = 4.2 * ms 
syns.append(rc_inh)
netcons.append(rc_inh_con)


