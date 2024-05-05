### Implementation of two models which study how positional information is encoded and interpreted in the developing embryo.

The models were developed alongside experiments performed in chick embryos, before the formation of the primitive streak. The site of the initiation of streak formation breaks radial symmetry, defining the posterior of the embryo. The experiments involved grafting beads soaked in compounds which either induce or inhibit streak formation.

We model a ring of cells around the circumference of the embryo, in a region called the marginal zone which has been shown to be key for the formation of the primitive streak. Each cell has a defined concentration of streak-inducer and -inhibitor. Cells use these concentrations to make the binary decision to initiate the formation of a streak, or not.

Both models involve cells balancing their concentrations of inducer and inhibitor. The key differences for models A and B are that
<ol type="A">
  <li>cells assess their values of inducer/inhibitor autonomously without reference to their neighbours,</li>
  <li>cells compare their own values of inducer/inhibitor with those of their neighbours.</li>
</ol>

## Quickstart

Code has been tested with Python 3.12.2.
#### Download code
```
git clone https://github.com/catohaste/neighbourhood-streak.git
cd neighoubourhood-streak
```

#### Create and activate virtual environments with all required packages
```
python3 -m venv env
source env/bin/activate
python -m pip install -r requirements.txt
```

#### Run files
```
python main_signal_slope.py
```

#### Jupyter demo
```
python -m ipykernel install --user --name=env
jupyter notebook
```

Open 'demo.ipynb' and change kernel to the virtualenv created, if not already done so.

### Clean up
```
jupyter kernelspec uninstall env
deactivate
rm -r env
```
