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
  printf("USAGE: %s <rates> <outfile>\n",program_name());
endfunction

## check command line arguments
if (nargin<2)
  printf("ERROR: please provide rates and output file.\n");
  help();
  exit(1);
endif
ratesfilename = arg_list{1};
outfilename = arg_list{2};


############################################################
## parameters
##
starttime  = 1;         # start time
endtime    = 1e12;      # end time
tinc       = 1.02;       # time increment
startstate = 2;         # state with initial probability 1


############################################################
##
## functions

function fprintvec(fout,v,dim)
  for i=1:dim
    fprintf(fout," %e",v(i,1));
  endfor
endfunction

function printvec(v,dim)
  fprintvec(stdout,v,dim);
endfunction



############################################################
##
## main

printf("Load rates from file %s\n",ratesfilename);
R = load("-ascii",ratesfilename);
dim=size(R,1);

if (size(R,2)!=dim)
  printf("ERROR: Rates file has to contain a square matrix.\n");
  exit(1);
endif  

format "short" "e";


pi0 = zeros(dim,1);
pi0(startstate,1)=1;
printf("PI_0: ");
printvec(pi0,dim)
printf("\n");

## recalculate diagonal of R (as -rowsum)
R = R - diag(diag(R)); # set diagonal to 0
rowsums = R * ones(dim,1);
R = R - diag(rowsums);
R = transpose(R);

pie = expm(endtime*R)*pi0;
printf("PI_e: ");
printvec(pie,dim)
printf("\n");

#printf("Rates:\n");
#disp(R);

printf("\n");
printf("Compute distributions at times %e..%e (until convergence)\n",starttime,endtime);

fout = fopen (outfilename, "w");

time=starttime;
while(time<endtime)
  pi = expm(time*R)*pi0;
  fprintf(fout,"%e",time);
  fprintvec(fout,pi,dim)
  fprintf(fout,"\n");
  
  diff=transpose(abs(pi-pie))*ones(dim,1);
  # disp(diff)
  if (diff(1,1) < 1e-8)
    printf("Convergence at time %e.\n",time);
    break;
  endif

  time *= tinc;
endwhile

fclose(fout);
