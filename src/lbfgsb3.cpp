//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Rmath.h>
#include <Rcpp.h>
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
using namespace Rcpp;

extern "C" void setulb_(int *n, int *m, double *x, double *l, double *u,
			int *nbd, double *f, double *g, double *factr, double *pgtol,
			double *wa, int *iwa, int *itask, int *iprint,
			int *icsave, int *lsave, int *isave, double *dsave);

typedef double optimfn(int n, double *par, void *ex);

typedef void optimgr(int n, double *par, double *gr, void *ex);

extern "C" void lbfgsb3C_(int n, int lmm, double *x, double *lower,
			  double *upper, int *nbd, double *Fmin, optimfn fn,
			  optimgr gr, int *fail, void *ex, double factr,
			  double pgtol, int *fncount, int *grcount,
			  int maxit, char *msg, int trace, int iprint,
			  double xtol, double *g){
  // Optim compatible interface
  int itask= 2;
  // *Fmin=;
  double *lastx = new double[n];
  std::copy(&x[0],&x[0]+n,&lastx[0]);
  int nwa = 2*lmm*n + 11*lmm*lmm + 5*n + 8*lmm;
  double *wa= new double[nwa];
  int niwa = 3*n;
  int *iwa= new int[niwa];
  int icsave = 0;
  int lsave[4] = {0};
  int isave[44] = {0};
  int i=0;
  double dsave[29]= {0};
  // Initial setup
  int doExit=0;
  fncount[0]=0;
  grcount[0]=0;
  int itask2=0;
  double tmp;
  while (true){
    setulb_(&n, &lmm, x, lower, upper, nbd, Fmin, g, &factr, &pgtol,
	  wa, iwa, &itask, &iprint, &icsave, lsave, isave, dsave);
    switch (itask){
    case 4:
    case 20:
    case 21:
      // Calculate f and g
      Fmin[0] = fn(n, x, ex);
      fncount[0]++;
      gr(n, x, g, ex);
      grcount[0]++;
      break;
    case 1:
      // New x;
      if (maxit <= fncount[0]){
	itask2=28;
	itask=3; // Stop -- gives the right results and restores gradients
      } else {
	tmp = fabs(lastx[n-1]-x[n-1]);
	for (i=n-1;i--;){
	  tmp = max2(tmp,fabs(lastx[i]-x[i]));
	}
	if (tmp < xtol){
	  itask2=27;
	  itask=3; // Stop -- gives the right results and restores gradients
	}
      }
      break;
    default:
      doExit=1;
    }
    if (doExit) break;
  }
  if (itask2){
    itask=itask2;
  }
  // info <- list(task = task, itask = itask, lsave = lsave,
  //      icsave = icsave, dsave = dsave, isave = isave)
  fail[0]= itask;
  delete[] wa;
  delete[] iwa;
  delete[] lastx;
}

Environment grho;

List ev;

double gfn(int n, double *x, void *ex){
  Rcpp::NumericVector par(n);
  std::copy(&x[0], &x[0]+n, &par[0]);
  Function fn = as<Function>(ev["fn"]);
  double ret = as<double>(fn(par, grho));
  return ret;
}

void ggr(int n, double *x, double *gr, void *ex){
  Rcpp::NumericVector par(n), ret(n);
  std::copy(&x[0], &x[0]+n, &par[0]);
  Function grad = as<Function>(ev["gr"]);
  ret = grad(par, grho);
  std::copy(&ret[0], &ret[0]+n, &gr[0]);
}

//[[Rcpp::export]]
Rcpp::List lbfgsb3cpp(NumericVector par, Function fn, Function gr, NumericVector lower, NumericVector upper, List ctrl, Environment rho){
  Rcpp::List ret;
  ev["fn"] = fn;
  ev["gr"] = gr;
  Rcpp::NumericVector g(par.size());
  CharacterVector taskList(28);
  taskList[0]="NEW_X";
  taskList[1]="START";
  taskList[2]="STOP";
  taskList[3]="FG";//,  // 1-4
  taskList[4]="ABNORMAL_TERMINATION_IN_LNSRCH";
  taskList[5]="CONVERGENCE"; //5-6
  taskList[6]="CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL";//7
  taskList[7]="CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH";//8
  taskList[8]="ERROR: FTOL .LT. ZERO"; //9
  taskList[9]="ERROR: GTOL .LT. ZERO";//10
  taskList[10]="ERROR: INITIAL G .GE. ZERO"; //11
  taskList[11]="ERROR: INVALID NBD"; // 12
  taskList[12]="ERROR: N .LE. 0"; // 13
  taskList[13]="ERROR: NO FEASIBLE SOLUTION"; // 14
  taskList[14]="ERROR: STP .GT. STPMAX"; // 15
  taskList[15]="ERROR: STP .LT. STPMIN"; // 16
  taskList[16]="ERROR: STPMAX .LT. STPMIN"; // 17
  taskList[17]="ERROR: STPMIN .LT. ZERO"; // 18
  taskList[18]="ERROR: XTOL .LT. ZERO"; // 19
  taskList[19]="FG_LNSRCH"; // 20
  taskList[20]="FG_START"; // 21
  taskList[21]="RESTART_FROM_LNSRCH"; // 22
  taskList[22]="WARNING: ROUNDING ERRORS PREVENT PROGRESS"; // 23
  taskList[23]="WARNING: STP .eq. STPMAX"; // 24
  taskList[24]="WARNING: STP .eq. STPMIN"; // 25
  taskList[25]="WARNING: XTOL TEST SATISFIED"; //
  taskList[26] = "CONVERGENCE: Parameters differences below xtol";
  taskList[27] = "Maximum number of iterations reached";
  // CONV in 6, 7, 8; ERROR in 9-19; WARN in 23-26
  int trace = as<int>(ctrl["trace"]);
  double factr = as<double>(ctrl["factr"]);
  double pgtol = as<double>(ctrl["pgtol"]);
  double xtol = as<double>(ctrl["xtol"]);
  int lmm = as<int>(ctrl["lmm"]);
  int n = par.size();
  int maxit = as<int>(ctrl["maxit"]);
  int iprint=as<int>(ctrl["iprint"]);
  // double *g = new double[par.size()];
  double *low = new double[par.size()];
  if (lower.size() == 1){
    std::fill_n(&low[0],par.size(),lower[0]);
  } else if (lower.size() == par.size()){
    std::copy(lower.begin(),lower.end(),&low[0]);
  } else {
    delete [] low;
    stop("Lower bound must match the size of par or only have one element.");
  }
  double *up = new double[par.size()];
  if (upper.size() == 1){
    std::fill_n(&up[0],par.size(),upper[0]);
  } else if (upper.size() == par.size()){
    std::copy(upper.begin(),upper.end(),&up[0]);
  } else {
    delete [] low;
    delete [] up;
    stop("Upper bound must match the size of par or only have one element.");
  }
  double *x = new double[par.size()];
  std::copy(par.begin(),par.end(),&x[0]);
  int *nbd = new int[par.size()];
  int i;
  for (i = par.size();i--;){
    /*
	   nbd(i)=0 if x(i) is unbounded,
		  1 if x(i) has only a lower bound,
		  2 if x(i) has both lower and upper bounds,
		  3 if x(i) has only an upper bound.
    */    
    nbd[i] = 0;
    if (R_FINITE(low[i])) nbd[i] = 1;
    if (R_FINITE(up[i]))  nbd[i] = 3 - nbd[i];
  }
  double fmin=std::numeric_limits<double>::max();
  int fail = 0, fncount=0, grcount=0;
  grho=rho;
  //void *ex = (void*)rho; //Should work but use global instead.
  void *ex =NULL;
  char msg[120];
  lbfgsb3C_(n, lmm, x, low, up, nbd, &fmin, gfn, ggr,
	    &fail, ex, factr, pgtol, &fncount,
	    &grcount, maxit, msg, trace, iprint , xtol, &g[0]);
  NumericVector parf(par.size());
  std::copy(&x[0],&x[0]+par.size(),parf.begin());
  ret["par"]=parf;
  ret["grad"]=g;
  ret["value"] = fmin;
  IntegerVector cnt = IntegerVector::create(fncount,grcount);
  ret["counts"] = cnt;
  switch (fail){
  case 6:
  case 7:
  case 8:
  case 27:
    ret["convergence"]=0;
    break;
  case 28:
    ret["convergence"]=1;
    break;
  case 23:
  case 24:
  case 25:
  case 26:
    ret["convergence"] = 51;
    break;
  case 9:
  case 10:
  case 11:
  case 12:
  case 13:
  case 14:
  case 15:
  case 16:
  case 17:
  case 18:
  case 19:
    ret["convergence"] = 52;
    break;
  }
  ret["message"]= CharacterVector::create(taskList[fail-1]);
  delete [] x;
  delete [] low;
  delete [] up;
  delete [] nbd;
  return ret;
}
