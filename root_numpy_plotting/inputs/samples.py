import os
user = os.environ['USER']

#TOP_DIRECTORY = "/tmp/" + user + "/HmumuTuplesVer_2_4/"
TOP_DIRECTORY = "/eos/user/l/luadamek/EOPTuples/haddedTogether/EOPTuples_v05/"

INPUT = \
{ "PythiaJetJet": ("Pythia8 LowMu and Dijet", [\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361020.jetjet.newextrap_newextrap_hist",\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361021.jetjet.newextrap_newextrap_hist",\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361022.jetjet.newextrap_newextrap1_hist",\
]),\
"LowMuData": ("2017 Data, Low-<#mu> Run 341294",[\
TOP_DIRECTORY + "user.luadamek.data17_13TeV.00341294.physics_MinBias.newextrap_newextrap_hist",\
])}
