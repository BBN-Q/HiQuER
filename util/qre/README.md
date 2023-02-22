# qre - quantum resource estimator

### Pre-requisites
- Python 3
- staq (https://github.com/softwareqinc/staq)

### Usage:

``` bash
cd src
cat file.qasm | staq_lattice_surgery | ./staq_qre.py -c config_compact.json
```

Or, for a given JSON file already produced by `staq_lattice_surgery`:

``` bash
cd src
cat file.json | ./staq_qre.py -c config_fast.json
```

To read the JSON directly (without stdin redirection):
``` bash
cd src
./staq_qre.py -c config_fast.json -f file.json
```

For help
``` bash
cd src
./staq_qre.py -h
```
