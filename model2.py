# =============================================================================
# MODEL 2 (3 RGC INPUTS AND TRIADS)
# -----------------------------------------------------------------------------
# Patrick McCarthy, Laura Lazzari and Jonathan Martin, March 2020 
# ----------------------------------------------------------------------------- 
# This file is to be used in conjunction with 'ballandsticks2.py'. It 
# implements classes from this file to create cell objects and then adds 
# synapses to simulate both the input stimuli and excitatory and inhibitory 
# connections between cells.
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

# stimulator 1
stim1 = h.NetStim() 
stim1.number = 1
stim1.start = 5 * ms 
stim1.interval = 0.5 * ms

# stimulator 2
stim2 = h.NetStim() 
stim2.number = 1
stim2.start = 5 * ms 
stim2.interval = 0.5 * ms

# stimulator 2
stim3 = h.NetStim() 
stim3.number = 1
stim3.start = 5 * ms 
stim3.interval = 0.5 * ms

# SYNAPSES

# set 1

# excitatory synapse on relay cell - technically part of triad but named differently for clarity
rc_exc1 = h.Exp2Syn(relaycell.soma(0.28))
rc_exc1_con = h.NetCon(stim1, rc_exc1) 
rc_exc1_con.delay = 0 * ms
rc_exc1_con.weight[0] = 5 # Hieberg: 11.6
rc_exc1.e = 42 * mV
rc_exc1.tau1 = 1 * ms 
rc_exc1.tau2 = 2 * ms 
syns.append(rc_exc1)
netcons.append(rc_exc1_con)

# proximal excitatory axodendritic synapse on interneuron
in_exc1 = h.Exp2Syn(interneuron.dend1_p(0.1)) 
in_exc1_con = h.NetCon(stim1, in_exc1) 
in_exc1_con.weight[0] = 0.6
in_exc1_con.delay = 0 * ms 
in_exc1.e = 42 * mV
in_exc1.tau1 = 1.6 * ms 
in_exc1.tau2 = 3.6* ms 
syns.append(in_exc1)
netcons.append(in_exc1_con)

# distal excitatory axodendritic synapse on interneuron 
triad1_1 = h.Exp2Syn(interneuron.dend1_d(1)) 
triad1_1_con = h.NetCon(stim1, triad1_1) 
triad1_1_con.weight[0] = 2
triad1_1_con.delay = 0 * ms 
triad1_1.e = 42 * mV
triad1_1.tau1 = 1 * ms 
triad1_1.tau2 = 2 * ms 
syns.append(triad1_1)
netcons.append(triad1_1_con)

# inhibitory dendrodendritic synapse on relay cell
triad1_2 = h.Exp2Syn(relaycell.soma(0.3)) 
triad1_2_con = h.NetCon(interneuron.dend1_d(0.99)._ref_v, triad1_2, sec=interneuron.dend1_d) 
triad1_2_con.weight[0] = 10
triad1_2_con.delay = 0.5 * ms 
triad1_2.e = -75 * mV
triad1_2.tau2 = 0.7 * ms 
triad1_2.tau2 = 4.2 * ms 
syns.append(triad1_2)
netcons.append(triad1_2_con)

# set 2

# excitatory synapse on relay cell - technically part of triad but named differently for clarity
rc_exc2 = h.Exp2Syn(relaycell.soma(0.58))
rc_exc2_con = h.NetCon(stim2, rc_exc2) 
rc_exc2_con.delay = 0 * ms
rc_exc2_con.weight[0] = 5 # Hieberg: 11.6
rc_exc2.e = 42 * mV
rc_exc2.tau1 = 1 * ms 
rc_exc2.tau2 = 2 * ms 
syns.append(rc_exc2)
netcons.append(rc_exc2_con)

# proximal excitatory axodendritic synapse on interneuron
in_exc2 = h.Exp2Syn(interneuron.dend2_p(0.1)) 
in_exc2_con = h.NetCon(stim2, in_exc2) 
in_exc2_con.weight[0] = 0.6
in_exc2_con.delay = 0 * ms 
in_exc2.e = 42 * mV
in_exc2.tau1 = 1.6 * ms 
in_exc2.tau2 = 3.6* ms 
syns.append(in_exc2)
netcons.append(in_exc2_con)

# distal excitatory axodendritic synapse on interneuron 
triad2_1 = h.Exp2Syn(interneuron.dend2_d(1)) 
triad2_1_con = h.NetCon(stim2, triad2_1) 
triad2_1_con.weight[0] = 2
triad2_1_con.delay = 0 * ms 
triad2_1.e = 42 * mV
triad2_1.tau1 = 1 * ms 
triad2_1.tau2 = 2 * ms 
syns.append(triad2_1)
netcons.append(triad2_1_con)

# inhibitory dendrodendritic synapse on relay cell 
triad2_2 = h.Exp2Syn(relaycell.soma(0.6)) 
triad2_2_con = h.NetCon(interneuron.dend2_d(0.99)._ref_v, triad2_2, sec=interneuron.dend2_d) 
triad2_2_con.weight[0] = 10
triad2_2_con.delay = 0.5 * ms 
triad2_2.e = -75 * mV
triad2_2.tau2 = 0.7 * ms 
triad2_2.tau2 = 4.2 * ms 
syns.append(triad2_2)
netcons.append(triad2_2_con)

# set 3

# excitatory synapse on relay cell - technically part of triad but named differently for clarity
rc_exc3 = h.Exp2Syn(relaycell.soma(0.88))
rc_exc3_con = h.NetCon(stim3, rc_exc3) 
rc_exc3_con.delay = 0 * ms
rc_exc3_con.weight[0] = 5 # Hieberg: 11.6
rc_exc3.e = 42 * mV
rc_exc3.tau1 = 1 * ms 
rc_exc3.tau2 = 2 * ms 
syns.append(rc_exc3)
netcons.append(rc_exc3_con)

# proximal excitatory axodendritic synapse on interneuron
in_exc3 = h.Exp2Syn(interneuron.dend3_p(0.1)) 
in_exc3_con = h.NetCon(stim3, in_exc3) 
in_exc3_con.weight[0] = 0.6
in_exc3_con.delay = 0 * ms 
in_exc3.e = 42 * mV
in_exc3.tau1 = 1.6 * ms 
in_exc3.tau2 = 3.6* ms 
syns.append(in_exc3)
netcons.append(in_exc3_con)

# distal excitatory axodendritic synapse on interneuron 
triad3_1 = h.Exp2Syn(interneuron.dend3_d(1)) 
triad3_1_con = h.NetCon(stim3, triad3_1) 
triad3_1_con.weight[0] = 2
triad3_1_con.delay = 0 * ms 
triad3_1.e = 42 * mV
triad3_1.tau1 = 1 * ms 
triad3_1.tau2 = 2 * ms 
syns.append(triad3_1)
netcons.append(triad3_1_con)

# inhibitory dendrodendritic synapse on relay cell 
triad3_2 = h.Exp2Syn(relaycell.soma(0.9)) 
triad3_2_con = h.NetCon(interneuron.dend3_d(0.99)._ref_v, triad3_2, sec=interneuron.dend3_d) 
triad3_2_con.weight[0] = 10
triad3_2_con.delay = 0.5 * ms 
triad3_2.e = -75 * mV
triad3_2.tau2 = 0.7 * ms 
triad3_2.tau2 = 4.2 * ms 
syns.append(triad3_2)
netcons.append(triad3_2_con)

# inhibitory axosomatic synapse between interneuron and relay cell
rc_inh = h.Exp2Syn(relaycell.soma(0.1)) 
rc_inh_con = h.NetCon(interneuron.axon_d(1)._ref_v, rc_inh, sec=interneuron.axon_d) 
rc_inh_con.weight[0] = 10
rc_inh_con.delay = 1  * ms 
rc_inh.e = -75 * mV
rc_inh.tau1 = 0.7 * ms 
rc_inh.tau2 = 4.2 * ms 
syns.append(rc_inh)
netcons.append(rc_inh_con)


