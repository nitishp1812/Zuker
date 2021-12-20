# Zuker
CS 466 Final Project

## Installation and Running

You need to have python3 installed to run the benchmarking file and the zuker algorithm

Install the required files using the following command:
```bash
pip3 install -r requirements.txt
```

### Zuker algorithm
You can run the Zuker algorithm using Python3:
```
python3 zuker.py --seq=GCGGAUUCUGCCCCAAUUCGCACCA
```

The seq command can be modified to get the output for any sequence as required.

### Benchmarking
To run the benchmarking file, you will need to download the [dot bracket files from the bpRNA-1m dataset](http://bprna.cgrb.oregonstate.edu/download.php#bpRNA-1m)
The zipped folder needs to be stored in the same directory as this repository

After this, you can run the benchmarking file using Python3
```
python3 benchmark.py
```

This file generates a file called output.txt and appends new lines to it in following runs. Delete the file if you would like a clean file to be created.

Currently, the benchmark files run a random sample containing approximately 10% of the data (which takes 1-1.5 hours to run). You can edit this percentage as needed on line 26 in `benchmark.py`.
