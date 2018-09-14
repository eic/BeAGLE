*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see LPTOU2)
*
      COMMON /py6strf/  py6f1, py6f2, py6R, sigobs, sigtrue
      DOUBLE PRECISION py6f1, py6f2, py6R, sigobs, sigtrue
      SAVE /py6strf/
