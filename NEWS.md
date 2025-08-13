## survivl 0.4.1
-------------------------------------------------------------------------------

CHANGES

 * Created `survivl_model`, a mirror image of the `causl_model` object.

 * New method inspired by the empirical rank technique from `causl` package and 
 the linear regression to residualize away correlations.
 
 * Switched all `msm_samp` to `survivl_model` and then calling `rmsm(n, surv_model)`

 * Switched some ipw calculations from by hand to using `ipw::ipwtm`.
 
 * Using `coxph` instead of `svyglm`.



## survivl 0.3.2
-------------------------------------------------------------------------------

CHANGES

 * Minor edits to facilitate censoring.



## survivl 0.3.1
-------------------------------------------------------------------------------

CHANGES

 * Introduced "weibull" survival times.

 

## survivl 0.3.0
-------------------------------------------------------------------------------

CHANGES

 * Substantial reorganization of `process_input()` code, partly to make more 
 use of `causl` package, and partly to  simplify by delegating tasks to other 
 functions.

 * Many function names changed to use `snake_case` rather than `camelCase`.
 
 
BUG FIXES

 * Fixed bug in `msm_samp` for multiple outcomes, so that survival was not 
 correctly established.
 
