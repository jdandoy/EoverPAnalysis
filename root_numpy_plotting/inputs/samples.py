import os
user = os.environ['USER']

#TOP_DIRECTORY = "/tmp/" + user + "/HmumuTuplesVer_2_4/"
TOP_DIRECTORY = "/Users/lukasadamek/cernbox/EOPTuples/haddedTogether/EOPTuples_v02/"

INPUT = \
{ "PythiaJetJet": ("Pythia Dijet", [\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361020.jetjet.Sept21_2018_EOP_noTrigger_hist",\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361021.jetjet.Sept21_2018_EOP_noTrigger_hist",\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361022.jetjet.Sept21_2018_EOP_noTrigger_hist",\
]),\
"MinBiasData": ("MinBias Data Run 00341294",[\
TOP_DIRECTORY + "user.luadamek.luadamek.data17_13TeV.00341294.Sept26_EOPTree_hist",\
])}
