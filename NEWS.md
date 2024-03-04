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
 
