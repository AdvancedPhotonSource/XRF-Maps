# XRF-Maps

## Running integrated spectra with custom preset

## 
```
--optimize-fit-override-params 0
```
Note that 0 means user specified custom preset from the override file

### Override file
For each fit parameter you can add the following 

fit_param_name + the following
```
fit_param_name + _MIN
fit_param_name + _MAX
fit_param_name + _STEPSIZE
fit_param_name + _FITTING
```

### _FITTING
```
fixed : do not fit this parameters
limit_lo_hi : If using mpfit, limit to the _MIN and _MAX ranges
limit_lo : If using mpfit, limit to the _MIN range
limit_hi : If using mpfit, limit to the _MAX range
fit : Fit this parameters with out an constraints.
```

### Example 
#### Energy offset fit it with limits
``` 
CAL_OFFSET_[E_OFFSET]:           -0.0054764490
CAL_OFFSET_[E_OFFSET]_MAX:        0.50000000
CAL_OFFSET_[E_OFFSET]_MIN:        -0.50000000
CAL_OFFSET_[E_OFFSET]_STEPSIZE:   0.001
CAL_OFFSET_[E_OFFSET]_FITTING:    limit_lo_hi
```

#### F_TAIL_OFFSET set limits but don't use them since we set FITTING to fit, will use step size if calling mpfit
```
F_TAIL_OFFSET:        0.003
F_TAIL_OFFSET_MAX:    0.100
F_TAIL_OFFSET_MIN:    0.00001
F_TAIL_OFFSET_STEPSIZE:    0.00001
F_TAIL_OFFSET_FITTING:    fit
```