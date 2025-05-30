## survivl 0.3.3.9002
-------------------------------------------------------------------------------

NEW FEATURES

 * Implemented joint family and link function processing.



## survivl 0.3.2.9000
-------------------------------------------------------------------------------

NEW FEATURES

 * Can now prespecify some data (beta feature)



## survivl 0.3.1.9000
-------------------------------------------------------------------------------

CHANGES

 * Introduced "weibull" survival times.

 

## survivl 0.3.0.9000
-------------------------------------------------------------------------------

CHANGES

 * Substantial reorganization of `process_input()` code, partly to make more 
 use of `causl` package, and partly to  simplify by delegating tasks to other 
 functions.

 * Many function names changed to use `snake_case` rather than `camelCase`.
 
 
BUG FIXES

 * Fixed bug in `msm_samp` for multiple outcomes, so that survival was not 
 correctly established.
 
