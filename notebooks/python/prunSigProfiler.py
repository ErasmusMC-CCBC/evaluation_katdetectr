from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerSimulator import SigProfilerSimulator as sigSim
from SigProfilerClusters import SigProfilerClusters as hp

sigSim.SigProfilerSimulator(
	project="results_SigProfiler",
	project_path="./",
	genome="GRCh37",
	contexts = ['288'],
	simulations=100,
	overlap=True)

hp.analysis(
	project="results_SigProfiler",
	genome="GRCh37",
	contexts="96",
	simContext=["288"],
	input_path="./",
	analysis="all",
	sortSims=True,
	subClassify=True,
	correction=True,
	calculateIMD=True,
	max_cpu=4,
	includedVAFs=False)
