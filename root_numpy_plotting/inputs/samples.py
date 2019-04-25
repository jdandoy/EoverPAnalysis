import os
user = os.environ['USER']

TOP_DIRECTORY="/eos/atlas/atlascerngroupdisk/perf-jets/EoverP/v01_tuples/"

INPUT = \
{ "PythiaJetJet": ("Pythia8X MinBias (JZ0W) and Dijet (361021 and JZ2W)", [\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist",\
]),\
"LowMuData": ("2017 Low-<#mu> Data",[\
TOP_DIRECTORY + "user.luadamek.data17_13TeV.00341294.physics_MinBias.calibhits_v01_hist",\
TOP_DIRECTORY + "user.luadamek.data17_13TeV.00341312.physics_MinBias.calibhits_v01_hist",\
TOP_DIRECTORY + "user.luadamek.data17_13TeV.00341419.physics_MinBias.calibhits_v01_hist",\
TOP_DIRECTORY + "user.luadamek.data17_13TeV.00341534.physics_MinBias.calibhits_v01_hist",\
TOP_DIRECTORY + "user.luadamek.data17_13TeV.00341615.physics_MinBias.calibhits_v01_hist",\
TOP_DIRECTORY + "user.luadamek.data17_13TeV.00341649.physics_MinBias.calibhits_v01_hist",\
])}
