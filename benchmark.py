import os
from tqdm import tqdm
import RNA
import random

RNA.params_load_RNA_Turner2004()

basepath = os.path.join('.', 'dbnFiles')
filenames = os.listdir(basepath)

datapoints = []

print('Loading datapoints...')
for filename in tqdm(filenames):
    filepath = os.path.join(basepath, filename)
    with open(filepath) as f:
        lines = f.readlines()

    sequence = lines[3].strip()
    dbn = lines[4].strip()
    
    datapoints.append((sequence, dbn, filename))

outFilepath = os.path.join('.', 'output.txt')

sample = random.sample(datapoints, int(0.1*len(datapoints)))
diff = 0.0
print('Calculating MFE and Zuker RNA folding for each datapoint')
for datapoint in tqdm(sample):
    seq, dbn, filename = datapoint

    zukerFold, zukerMFE = RNA.fold(seq)
    dataMFE = RNA.eval_structure_simple(seq, dbn)

    diff += zukerMFE - dataMFE

    if not os.path.exists(outFilepath):
        dataFormat = ['Format:\n', 'Sequence\n', 'Zuker (ViennaRNA) folding\n', 'Dataset folding\n', 'MFE for both structures\n\n\n']

        with open(outFilepath, 'w+') as f:
            f.writelines(dataFormat)

    with open(outFilepath, 'a') as f:
        data = [f'{seq}\n', f'{zukerFold}\n', f'{dbn}\n', f'Zuker MFE: {zukerMFE}\tDataset structure MFE: {dataMFE}\n\n\n']
        f.writelines(data)

avg = diff / len(sample)
print(f'Average difference: {avg} kcal/mol')
