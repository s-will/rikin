#!/usr/bin/octave -qf
# -*- Octave -*-

#
# INPUT: file with square matrix of partition functions of 1)
# transitions and 2) basins (given in diagonal).
#
# OUTPUT: probabilities of states at requested times according to
# Arrhenius rates form partitition functions; time is increased
# exponentially
#


tool_name="XRates";
tool_version="0.2";
tool_description="Calculate kinetics from weights/pfs of states and transitions by solving the master equation for Arrhenius rates.";

# ## like perls findbin, get the directory in that the script resides.
# ## e.g., this can be added to the function search path by addpath
# function ans = findbin()
#   ans = fileparts(make_absolute_filename(program_invocation_name));
# endfunction

# # printf("Script resides in directory %s\n",findbin());
# addpath(findbin());

############################################################
## simple octave option parser


function parser=argparse_addOption(parser,
			  option_name,
			  option_shortname,
			  option_type,
			  option_required,
			  option_help)
  
  n=length(parser);
  n++;

  if (length(option_name)==0)
    printf("argparse: non-empty (long) names required for all options.\n");
    return;
  endif

  positional=true;
  if ((length(option_name)>2 && substr(option_name,1,2)=="--"))
    positional=false;
    option_name=substr(option_name,3,length(option_name)-2);
  endif
  
  if length(option_shortname)>0
    positional = false;
  endif

  parser(n).positional=positional;
  parser(n).name=option_name;
  parser(n).shortname=option_shortname;
  parser(n).type=option_type;
  parser(n).required=option_required;
  parser(n).help=option_help;
endfunction

function p=argparse_parser
  p=[];
  p=argparse_addOption(p,"--help","h","",false,"Print this help.");
endfunction

function res=argparse_wrapText(text,prefix,width)
  splitat=0;
  splitskip=1;
  idx=strchr(text,"\n");
  if idx < width
    splitat=idx;
  elseif length(text)>width
    idx=strchr(substr(text,1,width)," ",1,"last");
    if (length(idx)>0)
      splitat=idx;
    else
      splitat=width;
      splitskip=0;
    endif
  endif
  
  if splitat>0
    leading_text=substr(text,1,splitat-splitskip);
    remaining_text = substr(text,splitat+1,length(text)-splitat);
    
    res = cstrcat(leading_text,
		  "\n",prefix,
		  strtrim(argparse_wrapText(remaining_text,prefix,width)));
  else
    res=text;
  endif
  
endfunction

function ans=argparse_buildOptionHelpString(p)
  namestr="";

  if (!p.positional)
    namestr=strcat(namestr,"--");
  endif
  namestr=strcat(namestr,p.name);
  
  if (length(p.shortname)>0)
    if (length(p.name)>0)
      namestr=strcat(namestr,",");
    endif
    namestr=strcat(namestr,"-",p.shortname);
  endif
  
  if (length(p.type)>0)
    namestr=strcat(namestr," <",p.type,">");
  endif
  
  helpstr=p.help;
  if (p.required)
    helpstr=strcat(helpstr," (mandatory)");
  endif
  
  prefix=sprintf("  %-30s","");
  helpstr=argparse_wrapText(helpstr,prefix,terminal_size()(2)-34);

  ans=sprintf("  %-30s%s",namestr,helpstr);
  
endfunction
  
function argparse_help(parser)
  printf("\n");
  
  for i=1:length(parser)
    p=parser(i);
    printf("%s\n",argparse_buildOptionHelpString(p));
  endfor
endfunction

function args = argparse_parse(parser,def_args)
endfunction

function ans=argparse_toolHeader(tool_name,tool_version,tool_description)
  ans=sprintf("%s %s -- %s",
	      tool_name,
	      tool_version,
	      argparse_wrapText(tool_description,"  ",74));
endfunction

function [args,arg_list,error_msg]=argparse_parseArgs(parser,def_args,arg_list)
  args=struct();
  error_msg="";
  
  positionals=[];
  
  for i=1:length(parser)
    p=parser(i);
    if p.positional
      positionals=[positionals,i];
    else
      
      idx=find(strcmp(arg_list,strcat("--",p.name)));
      
      if length(idx)==0 && length(p.shortname)>0
       	sidx=find(strcmp(arg_list,strcat("-",p.shortname)));
	if length(sidx)>0
	  if length(idx)>0
	    error_msg = strcat(error_msg,sprintf("Long and short option given for --%s.\n",p.name));
	  endif
	  idx=sidx;
	endif
      endif
      
      ## catch multiple options
      if length(idx)>1
	error_msg = strcat(error_msg,sprintf("Multiple options --%s disallowed.\n",p.name));
	return;
      endif
      
      remove_idxs=[];
      
      if strcmp(p.type,"")
	flag=true;
       	noidx=find(strcmp(arg_list,strcat("--no-",p.name)));
	if length(noidx)>0
	  if length(idx)>0
	    error_msg = strcat(error_msg,sprintf("Clash between arguments --%s and --no-%s.\n",p.name,p.name));
	  endif
	  remove_idxs=idx;
	  idx=noidx;
	  flag=false;
	endif
      endif

      if length(idx)!=0
	if (!strcmp(p.type,""))
	  
	  if length(arg_list)>idx
	    args.(p.name) = arg_list{idx+1};
	    
	    if (! strcmp(p.type,"string"))
	      args.(p.name) = str2num(args.(p.name));
	      if length(args.(p.name))==0
		error_msg = strcat(error_msg,sprintf("Value missing or of wrong type for option %s.\n",p.name));
		return;
	      endif
	    endif
	    
	    remove_idxs=[remove_idxs,idx+1];
	  else
	    error_msg = strcat(error_msg,sprintf("Value for option --%s missing.\n",p.name));
	  endif
	else 
	  args.(p.name) = flag;
	endif
	remove_idxs=[remove_idxs,idx];
      endif
      arg_list(remove_idxs,:)=[];
    endif
  endfor
  
  ## check for unknown options
  for i=1:length(arg_list)
    a=arg_list{i};
    if (length(a)>0 && strcmp(substr(a,1,1),"-"))
      error_msg=strcat(error_msg,sprintf("Unknown option %s.\n",a));
    endif
  endfor

  ## process positionals
  for i=1:length(positionals)
    p=parser(positionals(i));
    
    if (length(arg_list)==0) 
      if p.required
	error_msg = strcat(error_msg,sprintf("Positional argument '%s' is missing.\n",p.name));
      endif
    else
      args.(p.name) = arg_list{1};
      arg_list(1,:)=[];
    endif
  endfor

  ## check for required options
  for i=1:length(parser)
    p=parser(i);
    if p.required
      if !p.positional && ! isfield(args,p.name)
	error_msg = strcat(error_msg,sprintf("Required option '%s' is missing.\n",p.name));
      endif
    endif
  endfor
  
  ## set default arguments
  
  def_fieldnames=fieldnames(def_args);
  for i=1:length(def_fieldnames)
    fieldname=def_fieldnames{i};
    if (!isfield(args,fieldname)) 
      args.(fieldname)=def_args.(fieldname);
    endif
  endfor
  
endfunction

## end octave option parser
############################################################

############################################################
## default arguments
##
def_args = struct(
		  "verbose", false,
		  "binary", true,
		  "t0", 1e-1,
		  "t8", 1e16,
		  "tinc", 1.2,
		  "mode", "expm",
		  "p0", 2,
		  "absorb", -1,
		  "mfpts", false
		  );

parser=argparse_parser();
parser=argparse_addOption(parser,"--verbose","v","",false,"Be verbose.");
parser=argparse_addOption(parser,"--binary","b","",false,"Binary input");
parser=argparse_addOption(parser,"--nonbinary","","",false,"Non-binary input");
parser=argparse_addOption(parser,"--out","","string",true,
			  "Name of output file. The distributions at each time point until convergence are written as a table to this file.");
parser=argparse_addOption(parser,"--t0","","double",false,"Start time");
parser=argparse_addOption(parser,"--t8","","double",false,"End time");
parser=argparse_addOption(parser,"--tinc","","double",false,"Time increment");
parser=argparse_addOption(parser,"--mode","","string",false,"Mode (expm or diag)");
parser=argparse_addOption(parser,"--p0","","int",false,"State with initial probability 1");
parser=argparse_addOption(parser,"--absorb","","int",false,"Absorbing state");

parser=argparse_addOption(parser,"--mfpts","","",false,"Compute mean first passage times");

parser=argparse_addOption(parser,"in","","string",true,
			  "Input file. Matrix of transition weights/partition functions; the diagonal contains state weights.");

## note: optional positional arguments are allowed, but should be at the end (otherwise chaos prevails)
## parser=argparse_addOption(parser,"out2","","string",false,"Second output file");


printf("%s\n\n",argparse_toolHeader(tool_name,tool_version,tool_description));

[args,arg_list,error_msg] = argparse_parseArgs(parser,def_args,argv());

if (isfield(args,"help"))
  argparse_help(parser);
  exit(0);
endif

if (length(error_msg)>0)
  printf("Error: %s\n",error_msg);
  argparse_help(parser);
  exit(-1);
endif

############################################################
## parameters
##

verbose = isfield(args,"verbose") && args.verbose;   # verbose output

if verbose
  printf("Arguments:\n");
  disp(args);
endif


pffilename = args.in;
outfilename = args.out;


starttime  = args.t0;            # start time
endtime    = args.t8;            # end time
tinc       = args.tinc;          # time increment
startstate = args.p0;            # state with initial probability 1

#mode="diag";        # use diagonalization (fast, but errorprone)
mode=args.mode;        # use expm (slow, more stable)

if ! strcmp(mode,"expm") && ! strcmp(mode,"diag")
    printf("Error: unsupported mode %s.\n",mode);
    exit(-1);
endif

#binary  = isfield(args,"binary")  && args.binary;    # whether to read pfs in binary format
binary  = !( isfield(args,"nonbinary")  && args.nonbinary);    # whether to read pfs in binary format

absorb = args.absorb;

mfpts  = isfield(args,"mfpts")  && args.mfpts;

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
  m=(m+transpose(m))/2;
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

##############################
# Read and check pfs, preprocess

if (!binary)
  printf("Read text input from file %s\n",pffilename);
  pfs = load("-ascii",pffilename)
  dim=size(pfs,1);
else
  printf("Read binary input from file %s\n",pffilename);
  fh=fopen(pffilename);
  dim = fread(fh,1,"uint64");
  pfs = fread(fh,[dim,dim],"double");
  fclose(fh);
endif

##############################ä
## get basin pfs (from input pfs matrix)
#
basin_pfs = diag(pfs);


##############################
## symmetry check
## calculate rates from the original and the transposed pf matrix
## If the input is symmetric, then the two matrices are equal.
## Check the similarity of the matrices by norm(R1-R2,2).
## 
R1 = pfs ./ repmat(basin_pfs, 1,dim);
R1 = R1 - diag(diag(R1)); # set diagonal to 0

R2 = transpose(pfs) ./ repmat(basin_pfs, 1,dim);
R2 = R2 - diag(diag(R2)); # set diagonal to 0

dev_symm = norm( R1-R2,2);
if (dev_symm>0.001)
  printf("WARNING: questionable input symmetry of input pfs: norm_2(R(pfs)-R(t(pfs)))=%g (This value should be almost 0.)\n",dev_symm);
endif
## end symmetry check
##############################

##############################
## The input has to be symmetric; everything else is wrong input.
## Thus, in any case, calculate perfectly symmetric pfs.
## Symmetrize!
pfs = symmetrize(pfs);

if (size(pfs,2)!=dim)
  printf("ERROR: pf file has to contain a square matrix.\n");
  exit(1);
endif


##############################
## calculate R (from input pfs matrix)
#
R = pfs ./ repmat(basin_pfs, 1,dim);

## recalculate diagonal of R (as -rowsum)
R = R - diag(diag(R)); # set diagonal to 0
rowsums = R * ones(dim,1);
R = R - diag(rowsums);


format short e;

##############################
## initial distribution
#
pi0 = zeros(dim,1);
pi0(startstate,1)=1;
if (verbose)
  printf("PI_0: ");
  printvec(pi0)
  printf("\n");
endif

##############################
## stationary distribution
#
# pi8 = expm(endtime*R)*pi0; ## pi at endtime
pi8 = basin_pfs / (ones(1,dim)*basin_pfs); ## pi after infinite time

## handle absorbing state (optionally)
if (absorb>0)
  pi8=zeros(dim,1);
  pi8(absorb,1)=1;
  R(absorb,absorb)=0; ## set outflow of absorbing state to 0
  for i=1:dim
      R(absorb,i)=0;
  endfor
endif

if (verbose)
  printf("PI_8: ");
  printvec(pi8);
  printf("\n");
endif

#printf("R: \n");
#disp(R)


if (verbose)
  printf("\n");
  printf("Compute distributions at times %e..%e (until convergence)\n",starttime,endtime);
endif


## optionally, compute mean first passage times
if (mfpts)
  printf("Compute mean first passage times.\n");
  
  timescale=1; ## set "unit time": one step in the discrete MP with
  ## transition probability matrix P
  
  ## compute matrix of transition probabilities from rates:

  # P = eye(dim)+timescale*R; # approximation (like treekin): beginning of
  # Taylor expansion - this is usually good enough (since rates are very small)
  #
  P = expm(timescale * R); ## exact, P at unit time is the matrix exponential of timescale*R
  
  ## # Weather in the Land of Oz - example
  ## dim=3
  ## P = [2,1,1;2,0,2;1,1,2]/4 
  ## pi8 = P^100;
  ## pi8 = transpose(pi8(1,1:dim)) 
  
  if (absorb>0)
    # case: compute mfpt for the single state <absorb>, i.e. absorption times
    
    ## Q is the matrix P without the row and column of the absorbing state 
    Q = P;
    Q(:,absorb)=[]; # remove column absorb
    Q(absorb,:)=[]; # remove row absorb
    
    # fundamental matrix N
    
    N = inv(eye(dim-1)-Q);
        
    # row sums to get absorbtion times from each state i
    t = N * ones(dim-1,1);
    
    # insert absorbtion time 0 for absorb state
    t = [t(1:(absorb-1),1);0;t(absorb:(dim-1),1)];
    
    t=t*timescale;

    printf("Absorption times (to state %d): ",absorb);
    printvec(t);
    printf("\n");
    
  else
    # case: compute all mfpts
      
    # stationary matrix W
    W=repmat(transpose(pi8),dim,1);


    # compute fundamental matrix Z
    Z = inv( eye(dim) - P + W );
    # INV_TEST = sum( diag(Z*(eye(dim)-P+W)) ) / dim;


    # calculate mean first passage times
    # m_ij = (z_jj - z_ij) / w_j
    # m_ij is the mean first passage time from i to j
    M = ( repmat(transpose(diag(Z)),dim,1) - Z ) ./ W;
    
    M = M * timescale;
    
    disp(M);

    printf("MFPTs starting from state %d: ",startstate);
    printvec( M(startstate,:) )
    printf("\n");

  endif
endif


R=transpose(R); # we need R_ij = rate from j to i


## if we use diagonalization to exponentiate do some precomputation 
if (mode=="diag")
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
  #symmR=symmetrize(symmR); # force symmetrizing can make things worse  
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
  control3 = eigvecs * expdiag(eigvals) * eigvecs_inv;
  # disp(control3);
  
  dev_diag=norm(symmR - eigvecs * expdiag(eigvals) * eigvecs_inv,2);

  printf("Control diagonalization: %g (value should be almost 0)\n",dev_diag);

  
  ## precompute sub products
  pre_left  = sqrPI_ * eigvecs;
  pre_right = eigvecs_inv * _sqrPI * pi0;
endif


## open the output file to write distributions pi_t
fout = fopen (outfilename, "w");

step=0;
time=starttime;


while(time<endtime)
  if (mode=="diag")
    ## diagonalization mode
    pi = pre_left * expdiag(time*eigvals) * pre_right;
  else
    ## direct exponentiation mode
    pi = expm(time*R)*pi0;
  endif
  fprintf(fout,"%e",time);
  fprintvec(fout,pi)
  fprintf(fout,"\n");

  pi_sum=ones(1,dim)*pi;
  #if (pi_sum > 1.05)
  #  printf("Probability sum greater 1 at time %g (sum=%g).\n",time,pi_sum);
  #  break;
  #endif
  
  diff=transpose(abs(pi-pi8))*ones(dim,1);
  #disp(diff)
  if (diff(1,1) < 5e-3)
    printf("Convergence at time %g.\n",time);
    break;
  endif
  
  time *= tinc;
  step++;
endwhile

printf("Computed distributions at %d time points\n",step);

fclose(fout);
