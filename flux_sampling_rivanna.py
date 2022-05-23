from cobra.io import read_sbml_model
from cobra.sampling import sample

sbml_files = glob.glob(('*.sbml'))
bacteria = []
for file in sbmls_files:
	str(file)
	bac_name = file.split('.')
	bac_name.pop()
	bacteria+=bac_name

models = []
for file in glob.glob('*.sbml'):
	models.append(read_sbml_model(file, low_memory = False))

recons = {k:v for k, v in zip(bacteria, modesl)}

for bacteria, model in recons.items():
	flux = sample(model, 100)
	flux.to_csv(str(bacteria+'.csv'))