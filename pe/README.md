# Parameter estimation instructions

Requirements:
- teobresums (parameterized branch): https://bitbucket.org/teobresums/teobresums/branch/Dali-parametrized
- lalsuite (param teob branch): https://git.ligo.org/rossella.gamba/lalsuite-jacopo/-/tree/teobresums-dali-parameterized?ref_type=heads
- bilby (with param teob branch): 


## Setting up the environment
On LHO:
* source the conda environment under `/home/rossella.gamba/.conda/envs/lalsuite-param-teob`
* source the lalsuite installation under `/home/rossella.gamba/src/lalsuite-param-teob/_inst`:
```bash
source /home/rossella.gamba/src/lalsuite-param-teob/_inst/etc/lalsuiterc
```
* export the PYTHONPATH:
```bash
cd /home/rossella.gamba/src/lalsuite-param-teob/lalsimulation/python/lalsimulation
export PYTHONPATH=$PYTHONPATH:$PWD
```

## Parameters
The deviation parameters currently available are:
* `delta_a6c`
* `delta_cN3LO`
* `delta_abhf`
* `delta_Mbhf`
* `delta_Alm_mrg` 
* `delta_Omglm_mrg` 
* `delta_alphalm0` 
* `delta_omglm0`
They have to *all* be set in the prior file, either to a fixed value or by choosing an appropriate prior distribution.

Note that if one desires to sample the `delta_*lm` quantities, a dictionary must also be supplied to the waveform generator specifying for which mode (singular!) the deviation is considered.
TODO: allow for multiple mode deviations, find a way to implement things in a smarter way without adding all modes by hand in the generator function parameters...

## Example

An example config file and prior files are available here.
To run a PE, simply setup the environment as above, modify the `accounting-user` in the config file and then do
```bash
bilby_pipe config.ini
```
And follow the further instructions.
