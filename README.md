# SPT-3G + ACT DR6 + Planck NPIPE Lensing Likelihood

This repository contains likelihood software for the ACT + SPT CMB lensing analysis. If you use this software and/or the associated data, please cite the following papers:

- [SPT and ACT Collaborations (2025), arXiv: 2504.xxxxx](https://arxiv.org/abs/2504.xxxxx)
- [Ge, Millea, Camphuis, Daley, Huang, Omori, Quan et al SPT-3G Collaboration (2024), arXiv: 2411.06000](https://arxiv.org/abs/2411.06000)
- [Madhavacheril, Qu, Sherwin, MacCrann, Li et al ACT Collaboration (2023), arXiv:2304.05203](https://arxiv.org/abs/2304.05203)
- [Qu, Sherwin, Madhavacheril, Han, Crowley et al ACT Collaboration (2023), arXiv:2304.05202](https://arxiv.org/abs/2304.05202)
- [MacCrann, Sherwin, Qu, Namikawa, Madhavacheril et al ACT Collaboration (2024), arXiv: 2304.05196](https://arxiv.org/abs/2304.05196)

In addition, if you use the SPT+ACT+Planck lensing combination variant from the likelihood, please also cite:
- [Carron, Mirmelstein, Lewis (2022), arXiv:2206.07773, JCAP09(2022)039](https://arxiv.org/abs/2206.07773)

## Chains

Are we releasing chains? If so, where?


A pre-release version of the chains are available [here](). Please make sure to read the README file.

## Step 1: Install

First clone this repository, then install:

    pip install .

Tests can be run using 

    python -m unittest discover -s act_dr6_spt_lenslike/tests

## Step 2: download and unpack data

### ACT
Download the likelihood data tarball for ACT DR6 lensing from [NASA's LAMBDA archive](https://lambda.gsfc.nasa.gov/product/act/actadv_prod_table.html).

Extract the tarball into the `act_dr6_spt_lenslike/data/` directory in the cloned repository such the directory `v1.2` is directly inside it. Only then should you proceed with the next steps.

### SPT

Where?
Do we also need 2pt correction files?

### Covariances

Where?
    
## Step 3: use in Python codes

### Generic Python likelihood

```
import act_dr6_spt_lenslike as apslike

variant = 'act_baseline'
lens_only = False # use True if not combining with any primary CMB data
like_corrections = True # should be False if lens_only is True

# Do this once
data_dict = apslike.load_data(variant,lens_only=lens_only,like_corrections=like_corrections)
# This dict will now have entries like `data_binned_clkk` (binned data vector), `cov`
# (covariance matrix) and `binmat_act` (binning matrix to be applied to a theory
# curve starting at ell=0).

# Get cl_kk, cl_tt, cl_ee, cl_te, cl_bb predictions from your Boltzmann code.
# These are the CMB lensing convergence spectra (not potential or deflection)
# as well as the TT, EE, TE, BB CMB spectra (needed for likelihood corrections)
# in uK^2 units. All of these are C_ell (not D_ell), no ell or 2pi factors.
# Then call
lnlike=apslike.generic_lnlike(data_dict,ell_kk,cl_kk,ell_cmb,cl_tt,cl_ee,cl_te,cl_bb)
```

### Cobaya likelihood

Your Cobaya YAML or dictionary should have an entry of this form

```
likelihood:
    act_dr6_spt_lenslike.ACTDR6LensLike:
        lens_only: False
        stop_at_error: True
        lmax: 4000
        variant: act_baseline
```

No other parameters need to be set. (e.g. do not manually set `like_corrections` or `no_like_corrections` here).
An example is provided in `XXX.yaml`. If, however, you are combining with
the ACT DR4 CMB 2-point power spectrum likelihood, you should also set `no_actlike_cmb_corrections: True`
(in addition to `lens_only: True` as described below). You do not need to do this if you are combining
with Planck CMB 2-point power spectrum likelihoods.

### Important parameters

- `variant` should be
    - `act_baseline` for the ACT-only lensing power spectrum with the baseline multipole range
    - `act_extended` for the ACT-only lensing power spectrum with the extended multipole range (L<1250)
    - `actplanck_baseline` for the ACT+Planck lensing power spectrum with the baseline multipole range
    - `actplanck_extended` for the ACT+Planck lensing power spectrum with the extended multipole range (L<1250)
    - `spt_3g` for the SPT-only lensing power spectrum
    - `actspt3g_baseline` for the ACT+SPT lensing power spectrum
    - `actspt3g_extended` for the ACT+SPT lensing power spectrum with the extended multipole range (L<1250) for ACT
    - `actplanckspt3g_baseline` for the ACT+SPT+Planck lensing power spectrum with the baseline multipole range
    - `actplanckspt3g_extended` for the ACT+SPT+Planck lensing power spectrum with the extended multipole range (L<1250) for ACT
- `lens_only` should be
    - False when combining with any primary CMB measurement
    - True when not combining with any primary CMB measurement

### Recommended theory accuracy

For CAMB calls, we recommend the following (or higher accuracy):
- `lmax`: 4000
- `lens_margin`:1250
- `lens_potential_accuracy`: 4
- `AccuracyBoost`:1
- `lSampleBoost`:1
- `lAccuracyBoost`:1
- `halofit_version`:`mead2016`
