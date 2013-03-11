#!/usr/bin/octave -qf
printf ("XRates -- version 0.1 --");

arg_list = argv ();

## echo argument list
for i = 1:nargin
  printf (" %s", arg_list{i});
endfor
printf ("\n");

## print help message
function help()
  printf("USAGE: %s <pfs> <outfile>\n",program_name());
endfunction

## check command line arguments
if (nargin<2)
  printf("ERROR: please provide pf and output file.\n");
  help();
  exit(1);
endif
pffilename = arg_list{1};
outfilename = arg_list{2};


############################################################
## parameters
##
starttime  = 1e-1;         # start time
endtime    = 1e12;         # end time
tinc       = 1.02;         # time increment
startstate = 2;            # state with initial probability 1

use_diagonalization=true;

verbose=true;

############################################################
##
## functions

function fprintvec(fout,v)
  fprintf(fout," %g",v);
endfunction

function printvec(v)
  fprintvec(stdout,v);
endfunction

function [v,d] = eigsort (v,d)
  [dd,ix] = sort (-diag (d));
  for i=1:length(d)
    d(i,i) = - dd(i);
  end;
  v = v(:,ix)
endfunction

## force symmetrize (like done in treekin MxDiagonalize)
function m = symmetrize(m)
  n=length(m);
  for i=1:n
    for j=i+1:n
      m(i,j) = ( m(i,j)+m(j,i) ) / 2;
      m(j,i) = m(i,j);
    endfor
  endfor
endfunction

function res = expdiag(m)
  res=diag(exp(diag(m)));
endfunction


############################################################
##
## main

if (verbose)
  printf("Load pfs from file %s\n",pffilename);
endif
pfs = load("-ascii",pffilename);
# disp(pfs)

dim=size(pfs,1);
if (size(pfs,2)!=dim)
  printf("ERROR: pf file has to contain a square matrix.\n");
  exit(1);
endif

basin_pfs = diag(pfs);

R = pfs ./ repmat(basin_pfs, 1,dim);

## recalculate diagonal of R (as -rowsum)
R = R - diag(diag(R)); # set diagonal to 0
rowsums = R * ones(dim,1);
R = R - diag(rowsums);

#printf("R: \n");
#disp(R)

format "short";

pi0 = zeros(dim,1);
pi0(startstate,1)=1;
if (verbose)
  printf("PI_0: ");
  printvec(pi0)
  printf("\n");
endif

# pi8 = expm(endtime*R)*pi0; ## pi at endtime
pi8 = basin_pfs / (ones(1,dim)*basin_pfs); ## pi at endtime
if (verbose)
  printf("PI_8: ");
  printvec(pi8)
  printf("\n");
endif


if (verbose)
  printf("\n");
  printf("Compute distributions at times %e..%e (until convergence)\n",starttime,endtime);
endif


R=transpose(R); # we need R_ij = rate from j to i

if (use_diagonalization)
  ## make symmetric (works for rate matrices in detailed balance)
  ## ATTENTION: this requires to know the correct pi8
  ##
  
  sqrPI_ = diag(sqrt(pi8));
  _sqrPI = diag(sqrt(pi8).^(-1));
  symmR=R;
  symmR = symmR + diag(ones(1,dim));              # translate matrix
  symmR = _sqrPI*symmR*sqrPI_;

  # printf("symmetrized rates:\n");
  # disp(symmR);
  # symmR=symmetrize(symmR);
  # printf("force-symmetrized rates:\n");
  # disp(symmR);

  [eigvecs,eigvals]=eig(symmR);
  
  #[eigvecs,eigvals]=eigsort(eigvecs,eigvals);

  if (verbose)
    #printf("Eigenvalues:\n");
    #disp(transpose(diag(eigvals)));
  endif

  # printf("Eigenvectors:\n");
  # disp(eigvecs);
  
  eigvecs_inv  = transpose(eigvecs);
  #eigvecs_inv = inverse(eigvecs);
  
  # compensate for translation of matrix
  eigvals = eigvals - diag(ones(1,dim));

  # printf("Compute symmetrized rates from evecs and evals again:\n");
  # control1 = eigvecs * eigvals * eigvecs_inv;
  # disp(control1);

  # printf("Compute desymmetrized rates from evecs and evals again:\n");
  # disp(sqrPI_*(symmR-diag(ones(1,dim)))*_sqrPI);

  # printf("direct matrix exponentiation:\n");
  # control2 = expm(symmR-diag(ones(1,dim)));
  # disp(control2);

  # printf("from diagonalization:\n");
  # control3 = eigvecs * expdiag(eigvals) * eigvecs_inv;
  # disp(control3);
  
  ## precompute sub products
  pre_left  = sqrPI_ * eigvecs;
  pre_right = eigvecs_inv * _sqrPI * pi0;
endif


## open the output file to write distributions pi_t
fout = fopen (outfilename, "w");

time=starttime;


while(time<endtime)
  if (use_diagonalization)
    pi = pre_left * expdiag(time*eigvals) * pre_right;
  else
    pi = expm(time*R)*pi0;
  endif
  fprintf(fout,"%e",time);
  fprintvec(fout,pi)
  fprintf(fout,"\n");
  
  pi_sum=ones(1,dim)*pi;
  if (pi_sum > 1.05)
    printf("Probability sum greater 1 at time %g (sum=%g).\n",time,pi_sum);
    break;
  endif
  
  diff=transpose(abs(pi-pi8))*ones(dim,1);
  #disp(diff)
  if (diff(1,1) < 1e-3)
    printf("Convergence at time %g.\n",time);
    break;
  endif
  
  time *= tinc;
endwhile

fclose(fout);
