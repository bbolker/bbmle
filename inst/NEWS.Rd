\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{bbmle News}
\encoding{UTF-8}

\section{Changes in version 1.0.25}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{mle2.options()} function is now exported
    }
  }    
  \subsection{CRAN COMPATIBILITY}{
    \itemize{
      \item vignettes will build without suggested packages
    }
  }
}

\section{Changes in version 1.0.24 (2021-08-06)}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item new \code{namedrop_args} argument to \code{mle2} (don't drop
      argument names within sub-lists of data)
      \item allow user-specified covariance matrix for importance
      sampling
    }
  }
}
\section{Changes in version 1.0.23 (2020-02-03)}{
  \itemize{
    \item suppress some warnings in tests for cross-platform/CRAN compatibility
  }
}
\section{Changes in version 1.0.22 (2019-12-19)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item intercept-only parameters were handled wrong (GH #8)
      \item fix error-handling (class vs. inherits) for R 4.0 compatibility
    }
  }
  \subsection{NEW FEATURES}{
     \itemize{
       \item \code{pop_pred_samp} to compute population prediction
       samples, tweaking non-positive-definite covariance matrices and using
       importance samples if necessary/specified. \strong{This is
	 still VERY experimental/in flux!}
       \item added \code{dnorm_n} function for Normal densities
       with the standard deviation parameter profiled out
       \item add \code{slnorm} function for log-Normal densities
       (median and mean only for now)
       \item IC tables now have an \code{as.data.frame} method (from
       @iago-pssjd on GitHub)
     }
   }
   \subsection{USER-VISIBLE CHANGES}{
     \itemize{
       \item the evaluated
       \code{start} and \code{fixed} components of the call are now stored
       in the \code{call} component (Daniel B. Stouffer)
     } % itemize
   } % user-visible changes
 } % section 1.0.22
\section{Changes in version 1.0.20 (2017-10-30)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item fixed buglet: flipped profile plot axes, confint for negative values
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{summary.mle2} is now exported, for use by other
      packages
      \item small fixes to AIC(c) methods 
    }
  }
}

\section{Changes in version 1.0.19 (2017-04-08)}{
  \itemize{
    \item fixed bug: evaluate \code{call$method} so that
    profiling/updating works within a function environment
    \item make AICtab smarter about NA values
    \item fix BIC bug (infinite recursion)
    \item hessian computation uses gradient function if provided
    \item basic continuation method implemented for profiling (with
    stubs for smarter methods)
    \item mle2 stores its calling environment for more flexibility
    when re-evaluating, e.g. in profiling (could lead to occasional
    surprises, e.g. if saving a fitted mle2 object with large objects
    in its calling environment)
  }
}

\section{Changes in version 1.0.18 (2016-02-11)}{
  \itemize{
    \item update slice functionality; allow for explicit ranges
    \item CRAN updates (import from base packages)
  }
}
\section{Changes in version 1.0.17 (2014-01-01)}{
  \itemize{
    \item new warning if ~dnorm is called with sd implicitly ==1
    \item some internal tweaking to slice functionality
    \item updated to allow for \code{MuMIn} v 1.10.0 changes
  }
}

\section{Changes in version 1.0.16 (2014-01-01)}{
  \itemize{
    \item fix \code{mnames} behaviour
    \item slight vignette cleanup
  }
}
\section{Changes in version 1.0.15 (2013-11-20)}{
  \itemize{
    \item add \code{logLik} option to IC tables, more use of \code{nobs} methods
    \item minor improvements to \code{slice} functionality
  }
}

\section{Changes in version 1.0.14 (2013-08-24)}{
  \itemize{
    \item more CRAN tweaks
  }
}
\section{Changes in version 1.0.13 (2013-08-22)}{
  \itemize{
    \item .Rbuildignore tweaks for CRAN
  }
}

\section{Changes in version 1.0.12 (2013-08-22)}{
  \itemize{
    \item vignette tweaks
    \item add Depends: R>=3.0.0
  }
}
\section{Changes in version 1.0.11 (2013-08-19)}{
  \itemize{
    \item add .Rbuildignore for CRAN tests
  }
}
\section{Changes in version 1.0.10 (2013-08-18)}{
  \itemize{
    \item adapt to updated optimx
    \item tweaks for CRAN compliance
  }
}
\section{Changes in version 1.0.9 (2013-06-23)}{
  \itemize{
    \item switch from aod to aods3 in vignettes
    \item switch vignette to knitr
   }
}

\section{Changes in version 1.0.8 (2013-04-23)}{
  \itemize{
     \item tweaks to \code{print.ICtab()}
   }
}

\section{Changes in version 1.0.7 }{
  \itemize{
     \item warn on convergence failure
   }
}

\section{Changes in version 1.0.6 }{
  \itemize{
     \item fixed (fairly obscure) bug in assigning data environment to
  gradient function; replace 'gr' in call with appropriately evaluated version
   }
}

\section{Changes in version 1.0.5.3 (2012-09-05)}{
  \itemize{
     \item changed some cat() to message()
   }
}

\section{Changes in version 1.0.5.2 (2012-07-29)}{
  \itemize{
=======
    \item add .Rbuildignore for CRAN tests
  }
}
\section{Changes in version 1.0.10 (2013-08-18)}{
  \itemize{
    \item adapt to updated optimx
    \item tweaks for CRAN compliance
  }
}
\section{Changes in version 1.0.9 (2013-06-23)}{
  \itemize{
    \item switch from aod to aods3 in vignettes
    \item switch vignette to knitr
   }
}

\section{Changes in version 1.0.8 (2013-04-23)}{
  \itemize{
     \item tweaks to \code{print.ICtab()}
   }
}

\section{Changes in version 1.0.7 }{
  \itemize{
     \item warn on convergence failure
   }
}

\section{Changes in version 1.0.6 }{
  \itemize{
     \item fixed (fairly obscure) bug in assigning data environment to
  gradient function; replace 'gr' in call with appropriately evaluated version
   }
}

\section{Changes in version 1.0.5.3 (2012-09-05)}{
  \itemize{
     \item changed some cat() to message()
   }
}

\section{Changes in version 1.0.5.2 (2012-07-29)}{
  \itemize{
>>>>>>> .r116
     \item remove BIC definitions (now unnecessary/should be adequately
     defined in core R)
     \item add explicit Depends: on stats4
     \item note that development optimx (2012.05.24+) now fails on
     profiling when reduced model is 1-dimensional (i.e. for a
     2-parameter model)
   }
}

\section{Changes in version 1.0.5.1 (2012-07-17)}{
  \itemize{
     \item remove spurious .RData file; version bump
   }
}

\section{Changes in version 1.0.5 (2012-05-15)}{
  \itemize{
    \item wrapped eigen() call to prevent failure of eigenvalue ratio
    summary when hessian is bad
    \item fix bug: forgot to export summary method for mle2
    \item add \code{exclude.fixed} argument to \code{coef} method
    \item fix bug: single lower/upper parameter, or
    prof.lower/prof.upper parameter, not interpreted properly in
    profile/confint
    \item add and document \code{slice} methods: change from old
    (profile-like) behavior, old version is available as \code{sliceOld}
    \item DESCRIPTION/NAMESPACE fixes, move most Depends: to Imports:
    instead (except methods package)
  }
}

\section{Changes in version 1.0.4.2 (2012-02-25)}{
  \itemize{
    \item fix bug in gradient/vecpar/profiling interaction (Daniel
    Kaschek)
    \item improve (and document) uniroot method for confint: now
    respects box constraints
    \item fix issue where bobyqa (optimx) strips parameter names from result
  }
}
\section{Changes in version 1.0.4.1 (2012-01-27)}{
  \itemize{
    \item remove a bit of installed junk; vignette fix
  }
}
\section{Changes in version 1.0.4 (2012-01-02)}{
  \itemize{
    \item fix imports for latest version of R-devel
  }
}
\section{Changes in version 1.0.3}{
  \itemize{
    \item modified starting value code and documentation slightly:
    now allows per-parameter lists.  Updated docs.
    \item Fixed bug that would screw things up if 'i' were used
    as a parameter.
    \item Update vignette to conform to MuMIn changes
  }
}

\section{Changes in version 1.0.2 (2011-09-07)}{
  \itemize{
    \item fixed buglet in detection of no-intercept models when setting
    starting values
  }
}
\section{Changes in version 1.0.1 (2011-08-04)}{
  \itemize{
    \item Turned off some warnings when skip.hessian=TRUE (e.g. in profile)
    \item Calculate max grad and Hessian eigenval ratio, add to "details"
    \item (Should add accessor methods)
    \item Fixed NAMESPACE to export methods etc etc properly; removed
    universal exportPattern directive
  }
}

\section{Changes in version 1.0.0 (2011-06-17)}{
  \itemize{
    \item mainly just a version bump for CRAN
    \item added tag for variables in profile/confint warnings
    \item lots more technical information in vignette about profiling algorithm
}}

\section{Changes in version 0.9.9 (2011-05-14)}{
  \itemize{
    \item changed NEWS to NEWS.Rd
    \item fixed bug for "optimize" -- profile etc. now respect bounds
    \item eliminated warning about bounds names if all identical
    \item add "try_harder" flag to profiling (ignore flat spots, NAs ...)
}}
  
\section{Changes in version 0.9.8}{
  \itemize{
    \item gradient functions work better with fixed parameters, hence with profiling
    \item profile plot reverts to linear for non-monotonic profile
    \item added warning in confint for non-monotonic profile, or for non-monotonic spline fit to monotonic profile; revert from spline+linear to linear approximation in this case
    \item various documentation improvements
    \item optimx improvements
    \item require data= argument when using formula interface
    \item turn off hessian computation in profile
    \item allow use of MASS::ginv
}}

\section{Changes in version 0.9.7}{
  \itemize{
    \item bug fix in calc_mle2_function for no-intercept models
    (thanks to Colin Kremer)
    \item fixed optimx, added 'user' option
}}

\section{Changes in version 0.9.6}{
  \itemize{
    \item changed hessian calculation to use numDeriv code (causes tiny changes to hessian results that could matter in edge cases). Too lazy to provide a backward compatibility mode ...
    \item documented optimizer= choices in ?mle2
}}

\section{Changes in version 0.9.5.1}{
  \itemize{
    \item fixed bug in AICc (David Harris)
}}

\section{Changes in version 0.9.5}{
  \itemize{
    \item added NAMESPACE, various fixes to go with that
    \item beginnings of an RUnit testing framework
    \item tweaked vignette
    \item added prof.lower, prof.upper to profile()
    \item added "optimize" to list of allowed optimizers, some bug fixes
}}

\section{Changes in version 0.9.4.1}{
  \itemize{
    \item tweaked par() resetting in profile plots
}}

\section{Changes in version 0.9.4}{
  \itemize{
    \item more qAICc fixing
}}

\section{Changes in version 0.9.3 (2009-09-18)}{
  \itemize{
    \item tweaked handling of bounds: profile now succeeds
    on some 1D problems where it didn't before
    \item added deviance, residuals methods
    \item added newparams argument to predict, simulate; newdata argument to simulate
    \item added vignette (stub)
    \item added explicit params argument, to help sort out full parameter specifications when parameters is non-NULL
}}

\section{Changes in version 0.9.2 (2009-08-10)}{
  \itemize{
    \item fixed predict() for case with parameters
    \item added snorm
    \item changed ICtab defaults to weight=TRUE, base=FALSE, sort=TRUE
}}

\section{Changes in version 0.9.1}{
  \itemize{
\item added simulate method (formula interface only)
\item fix AICctab bug
\item remove spurious cat/print in profile
\item fix qAIC bug
}}

\section{Changes in version 0.9.0 (2008-08-26)}{
  \itemize{
  \item fix Tom Hobbs bug: named lower/upper/parscale/ndeps
  get rearranged properly, otherwise rearrange in order
  of "start" and issue a warning
  \item documentation tweak for S4 as.data.frame
  \item added sbeta to list of known distributions
  \item removed nlme requirement & auto-loading
 }}

 \section{Changes in version 0.8.9 (2008-08-04)}{
   \itemize{
     \item version bump, submit to CRAN
     \item added predict method
 }}

 \section{Changes in version 0.8.8 (2008-07-10)}{
   \itemize{
\item added flexibility for profile plotting (main, x labels etc.);
    added examples
\item added an instance of "namedrop" to fix naming problem
\item added tol.newmin to slice etc.
\item added check for numeric return from profile within confint
\item fixed bugs in profile plotting when profile is restricted
   to a subset of variables
\item added tests for par() to reset to original on exit
\item improved profile documentation
\item replicate std.err if specified in profile
\item add as.data.frame 
\item tweak tol.newmin (better fit found during profile) code
}}

\section{Changes in version 0.8.7 (2008-05-12)}{
  \itemize{
    \item version bump, moved to R-forge.
    \item reordered NEWS file (most recent first)
}}
\section{Changes in version 0.8.6.1 (2008-03-22)}{
  \itemize{
    \item tweaked stop-on-better-fit code
    \item fixed (?) qAIC(c) methods
}}

\section{Changes in version 0.8.6 (2008-03-26)}{
  \itemize{
    \item tweak/fix to ICtab documentation (thanks to Tom Hobbs)
    \item added qAIC(c) methods (not working yet!)
}}

\section{Changes in version 0.8.5.1}{
\itemize{
  \item oops. Fixed infelicity (bug?) in new  environment manipulation
  }}

\section{Changes in version 0.8.5}{
  \itemize{
    \item tweaked environment/data assignment to preserve
    original minuslogl environment better
}}

\section{Changes in version 0.8.4}{
  \itemize{
    \item changed  plot.profile.mle2 options (added onepage etc.,
    made plot.confstr=TRUE by default)
}}

\section{Changes in version 0.8.3}{
  \itemize{
    \item added warning about too-short lower/upper 
    \item added documentation
    }}

\section{Changes in version 0.8.2}{
  \itemize{
\item   fixed bug in AICctab
\item cosmetic change to printing -- save call.orig
\item moved ChangeLog to NEWS
}}

\section{Changes in version 0.8.1}{
  \itemize{fixed (?) environment bug
\item tried to use built-in relist, but failed: renamed relist
	  to "relist2" (try again later)
\item documented get.mnames (auxiliary function for ICtabs)
\item started to add gr (gradient) capability -- NOT TESTED
}}
\section{Changes in version 0.8}{
  \itemize{
    \item changed ICtab to allow either ICtab(x,y,z) or ICtab(list(x,y,z))
    (L <- list(...); if is.list(L[[1]]) && length(L)==1)
}}

\section{Changes in version 0.7.7}{
  \itemize{
    \item fix bug in profiling: all optim() methods EXCEPT L-BFGS-B.
    return the value of the objective function if given a function
    with no arguments/zero-length starting parameter vector
    (this is the situation with "profiling" a 1-D function).
    L-BFGS-B gives funky answers.  added a check for this case.
    (may need to check behavior for alternate optimizers (nlm etc))
    [this behavior triggered a "found better fit" error when profiling
     1D functions with L-BFGS-B]

     \item changed behavior when finding better fit during profiling
     to return new parameters
 }}


\section{Changes in version 0.7.6}{
  \itemize{
    \item tweak vignette
    \item fixed second major AICc bug (was fixed in mle2 method,
    but not in logLik method)
}}
     
\section{Changes in version 0.7.5}{
  \itemize{
    \item change "ll" to "LL" in examples for clarity
    \item tweaked anova reporting of models (wrap instead of truncating)
    \item added (undocumented) show.points option to profile plot
    to display actual locations of profile evaluation
\item tweaked profile to behave better when profiling variables
	with constraints (upper, lower)
\item moved vignette to inst/doc where it belongs
\item ICtab hack to protect against package:aod definition of AIC(logLik)
\item added submit stub
\item tweaked slice.mle2-class docs for consistency
\item fiddled with vignette
\item preliminary code to allow non-monotonic profiles
\item preliminary add nlm to list of optimizers (untested)
\item add aod, Hmisc, emdbook to VignetteDepends and Suggests:
}}

\section{Changes in version 0.7}{
  \itemize{
    \item better df extraction in ICtab
\item minor bug fix for AICc (allows AICc of nls objects)
\item handle models with -1 in formula better:
	starting values set "all equal"
\item made ANOVA formula line-length accessible
\item added skip.hessian and trace arguments to mle2
\item messed around with BIC definition -- attempt at consistency with nlme
\item added rudimentary support for nlminb, constrOptim
\item nlme now required for fdHess (which is required for
	nlminb since it doesn't compute a finite-diff
	Hessian)
    }}
    
    \section{Changes in version 0.6}{
      \itemize{
	\item add experimental formula interface
\item change all names from mle to mle2 to avoid confusion/conflicts
\item with stats4 version of mle
\item change internal structure of data evaluation
\item worked on vignette
\item added optimizer slot (stub)
}}

\section{Changes in version 0.5}{
  \itemize{
  \item fix AICc bug!  (was deviance+2*k*(k+1)/(n-k-1), not
  AIC+2*k*(k+1)/(n-k-1))
  }}

  \section{Changes in version 0.4}{
      \itemize{
    \item change AIC to AICc for corrections
    \item add AICtab for weights, delta, sort ... options
\item expose error messages occuring within profile()
\item uniroot tries harder to find a valid endpoint
\item truncate terms in anova.mle at 80 characters
}}

\section{Changes in version 0.3}{
  \itemize{
    \item enhanced anova method, works with print.anova
    \item tweaked namedrop() code -- ??
}}

\section{Changes in version 0.2}{
  \itemize{
    \item added parnames, parnames<-
    \item minor fix to allow "profiles" of 1-parameter models
    (skip fdHess call)
    \item minor change to print method for mle results
    \item tweaking "vecpar" (to allow parameter vectors in objective function)
    \item removed fdHess/nlme dependency
  }
}
