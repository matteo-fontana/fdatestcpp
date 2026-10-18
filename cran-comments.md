## Submission

This is a new package.

## Test environments

* Windows 11 x64, R 4.6.1 (local; Rtools45, g++ 14.3; ggplot2 4.0.3):
  `R CMD check --as-cran`.
* The same, with the package compiled with `-Wall -Wextra -pedantic` (no
  warnings from the package's own code), with `-fsanitize=undefined` and
  `-D_GLIBCXX_ASSERTIONS` (tests and examples clean), and with long double
  replaced by double (as on macOS arm64; tests pass).

## R CMD check results

0 errors | 0 warnings | 2 notes

* "New submission".
* "Found the following files/directories: ''NULL''" in the check directory.
  This appears on the local machine for any package, including an empty
  dummy package, and is not created by this package.

The incoming spell check may list author names (Abramowicz, Pini, Vantini)
and "resamples"; these are correct.

## Relationship with 'fdatest'

'fdatestcpp' re-implements in C++ the interval-wise and threshold-wise testing
procedures of the development version of 'fdatest'
(<https://github.com/permaverse/fdatest>, commit 67bb41c, licensed GPL >= 3),
with faster algorithms (the interval-wise adjustment costs O(B p^2) instead
of O(B p^3)) and additional resampling schemes for dependent data. With the
same inputs and random seed it reproduces that version's p-values exactly.
The CRAN release of 'fdatest' (2.1.1, GPL-2) contains the earlier interval
testing procedure only, not these functions, and no code was taken from it.
The authors of 'fdatest' (Pini, Stamm and Vantini) are listed as contributors
and copyright holders in Authors@R. Function and argument names follow
'fdatest' so that the functions can be used as drop-in replacements.

'fdatest' is in Suggests only: the tests can compare results with its
development version, but these comparisons are opt-in (environment variable
FDATESTCPP_COMPARE_FDATEST=true) and never run on CRAN.

Results have classes c("fdatestcpp_fos", "fos") and c("fdatestcpp_fts", "fts"):
"fos" and "fts" are the classes of the development version of 'fdatest', kept
so that its functions accept these results. 'fdatestcpp' registers its plot()
and autoplot() methods only for its own classes, so it never overrides
methods of another package.
