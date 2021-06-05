/* A2A3:  A Newey-West GMM test of A2 and A3.
*/

@------Name of the Gauss data set.---------------------------------@

dataset = "ewdset";
outwidth  132;

@--------Open GAUSS data set and define variable indices.----------@

    closeall f1;
    open f1=^dataset;

@-------Set the year for estimation and name of file for results. -@

     output file = idtest.out reset;

@-----------------Set the maximum number of iterations.------------@

maxiter = 299;

@-----------------Set the limit on the number of squeezes.---------@

    maxsqez=240;

@------------------Set the number of years.-----------------------@

    nyr = 4;

@-----------------Set the number of observations.------------------@

    nco = rowsf(f1)/nyr;



@---------------Initialize the counters.---------------------------@

   csave = zeros(nyr*16,2);

@-------Set the number of equations and parameters and name them. -@


nc = 2;
neq = 2;

let cnames = Ey2x Eyx2;

    @ ------- Now we enter the main body of the program. ----------- @

  cf = 0;
  do while cf <= 15;

      if cf==0;  "Deviations from means";
  elseif cf==1;  "With cash flow";
  elseif cf==2;  "With cash flow and d";
  elseif cf==3;  "With cash flow and b";
  elseif cf==4;  "With cash flow and s";
  elseif cf==5;  "With cash flow and d plus";
  elseif cf==6;  "With cash flow and b plus";
  elseif cf==7;  "With cash flow and s plus";
  elseif cf==8;  "With cash flow and db";
  elseif cf==9;  "With cash flow and ds";
  elseif cf==10; "With cash flow and bs";
  elseif cf==11; "With cash flow and db plus";
  elseif cf==12; "With cash flow and ds plus";
  elseif cf==13; "With cash flow and bs plus";
  elseif cf==14; "With cash flow and dbs";
  elseif cf==15; "With cash flow and dbs plus";
  endif;

  year = 1;
  do while year <= nyr;

   format 4,2;
   "Year = ";; year;
   print;


    @ --------- Read in the data and define the variables. ----------@


      dta=readr(f1,nco);

      n=rows(dta);

    y = dta[.,4];
    x = dta[.,5];
    z = dta[.,6];
    icept=ones(n,1);


      ww1    = (1-(dta[.,10].==0));
      ww2    = (dta[.,9].==9999);



      big1   = sortc(dta[.,7],1);
      one3b1 = big1[ceil(1*rows(big1)/3),1];
      two3b1 = big1[ceil(2*rows(big1)/3),1];
      big2   = sortc(dta[.,15],1);
      one3b2 = big2[ceil(1*rows(big2)/3),1];
      two3b2 = big2[ceil(2*rows(big2)/3),1];
      ww3 = dta[.,7] .< one3b1 .and dta[.,15] .< one3b2;

      wz1 = ww1.*z;
      wz2 = ww2.*z;
      wz3 = ww3.*z;



    if cf==0;       iz = icept;
    elseif cf==1;   iz = icept~z;

    elseif cf==2;   iz = icept~z~wz1;
    elseif cf==3;   iz = icept~z~wz2;
    elseif cf==4;   iz = icept~z~wz3;

    elseif cf==5;   iz = icept~z~wz1~ww1;
    elseif cf==6;   iz = icept~z~wz2~ww2;
    elseif cf==7;   iz = icept~z~wz3~ww3;

    elseif cf==8;   iz = icept~z~wz1~wz2;
    elseif cf==9;   iz = icept~z~wz1~wz3;
    elseif cf==10;  iz = icept~z~wz2~wz3;

    elseif cf==11;  iz = icept~z~wz1~wz2~ww1~ww2;
    elseif cf==12;  iz = icept~z~wz1~wz3~ww1~ww3;
    elseif cf==13;  iz = icept~z~wz2~wz3~ww2~ww3;

    elseif cf==14;  iz = icept~z~wz1~wz2~wz3;
    elseif cf==15;  iz = icept~z~wz1~wz2~wz3~ww1~ww2~ww3;
    endif;

    muy = inv( moment(iz,0) )* (iz'y);
    mux = inv( moment(iz,0) )* (iz'x);

    y_d = y - iz*muy;
    x_d = x - iz*mux;


    @ --------- Construct the moments for subsequent estimation. ----@
    y2  = y_d^2;
    x2  = x_d^2;
    yx  = x_d.*y_d;
    y2x = (y_d^2).*x_d;
    yx2 = (x_d^2).*y_d;

    Ey2x = meanc(y2x);
    Eyx2 = meanc(yx2);


    inEzz   =   invpd(moment(iz,0)/n);
    Ey2_z   =   meanc(y2   .* iz);
    Eyx_z   =   meanc(yx   .* iz);
    Ex2_z   =   meanc(x2   .* iz);
    zy_d    =   iz.*y_d;
    zx_d    =   iz.*x_d;

    @ ------- Initialize the parameters. --------------------------- @


    c = zeros(nc,1);

    c[1,1] = Ey2x;
    c[2,1] = Eyx2;

    @ --------- Here we enter the two-round estimation procedure. ---@
    rlim = 1;
    r=1;

    do while r <= rlim;

     @ --------The program creates the weighting matrix. ------------@
     clear w;

  w = optw(c);

     @ ------------- Start of iteration loop ------------------------@
iter = 1;
    tol=1e-9;                 @ Set the convergence criterion. @
    dc=1;                     @ Initialize the step length. @
do until abs(dc) < tol;
    clear g, f, sqz;

    f=deff(c);                @ The program jumps to the subroutine that
                                defines the f vector. @

    g=grad(c);                @ The program jumps to the subroutine that
                                computes analytic derivatives. The matrix
                                of partials of f with respect to the
                                parameters is called g. @

    obj = f'w*f;              @ This computes the value of the objective
                                function.@

    gwg= g'w*g;               @ This uses the GAUSS-NEWTON method
                                to compute full step dc. @
    gwf = g'w*f;


    onemore = 1;
    solveit:
    trap 1;
    dc = solpd(gwf,gwg);
    trap 0;
    if scalerr(dc) > 0;
      print "didn't invert! ";; "iter = ";; iter;; "round = ";; r;
      if onemore == 1;
        call chol(error(2));
        onemore = 0;
        goto solveit;
      else;
        print "whoops! ";; "iter = ";; iter;; "round = ";; r;
        end;
      endif;
    endif;


 if maxsqez > 0;              @ This jumps to the subroutines that
                                adjust the step length. @
   { c_new,sqz } = squeeze(c,dc);
 else;
   c_new=c - dc;
 endif;


 dc=c_new-c;                     @ Update variables for the next iteration. @
 c=c_new;
/*
output off;
format 1,0;
if iter == 1; cls; endif;
locate 2,1;
" ********** Computing GMM Iterations, Round " ;; r;; " *************** ";
?;

format 4,2;                   @ Print results for this iteration. @

" iter=";; iter;;
" No. of squeezes req'd=";; sqz;;
format 10,6;
" Starting value of obj=" nco*obj; @ This is the value of the objective function
                                   at the previous parameter estimates. @

@---------------Print out the results for the current iteration------------@

"          Coef.         Value         Step";
 format 14,6;
 mm=1;
 do until mm > nc;
  $cnames[mm,1];; "  ";; c[mm,1];; "  ";; dc[mm,1];
  mm=mm+1;
  endo;
*/
 iter=iter +1;
                                         @ Quit iterating if necessary. @
 if iter >= maxiter; goto escape; endif;
 if key == 27; goto escape; endif;       @ Hitting the escape key will stop
                                           the iteration loop.  @
endo;
@ ------ End of iteration loop ------------------------------------ @
escape:

@ Compute t-ratios, etc.                                      @

if r == rlim;
  w = optw(c);
  g=grad(c);
  gwg = g'w*g;
endif;

clear vc;
vc=invpd(gwg)/n;
stderr=sqrt(diag(vc));
t=c./stderr;
df=n-nc;
pvt=2*cdftc(abs(t),df);
@ -------------------------------------------------------- @
@ Print results.                                           @

cls;
if r==rlim;
output on;
endif;
/*
     "         GENERALIZED METHOD OF MOMENTS RESULTS, ROUND"  r;
     print;
     "Observations:"  nco;; "     Degrees of freedom:" df;
     print;
     "              Value of Objective Function: " nco*obj;
     print;
"         Coef       Value         Std. Error     T-Stat       P-Value";

     i=1;
     do until i > nc;
        format 14,6;
        $cnames[i,1];; c[i,1];; stderr[i,1];; t[i,1];; /rd pvt[i,1];
     i=i+1;
     endo;
*/
clear g, f, sqz, gwg, gwf, dc, stderr, t, pvt;

r=r+1;
endo;  @ This ends the "round" loop. @

     print;

nww = n*c'w*c;
print "Newey-West Wald statistic =" nww;
test=cdfchic(nww,2);
print;
print "P-value =" test;

csave[nyr*cf+year,1] = nww;
csave[nyr*cf+year,2] = test;

year = year + 1;
endo;

call seekr(f1,1);
cf = cf + 1;
endo;



/*--------------This business makes a LaTeX table------------*/

format /rd 8,3;


for qq(0, 15, 1);
  for jj(1, nyr, 1); "& ";; csave[nyr*qq+jj,1];; "  ";; endfor; " \\" "\\";
  for jj(1, nyr, 1); "&(";; csave[nyr*qq+jj,2];; ") ";; endfor; " \\" "\\";
endfor;


goto eop;  @ This skips all of the subroutines.        @

@ ----------------- Subroutines follow --------------- @

@ Subroutine to define the f vector. @

proc deff(a);
local f;

f = zeros(neq,1);

f[1,1] = Ey2x - a[1,1];
f[2,1] = Eyx2 - a[2,1];

retp(f);
clear f;
endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad(a);

retp(-eye(nc));

endp;

@  ------------------------------------------------------- @
@ Subroutine to compute squeezes.

  This subroutine compares the values of the objective function at the
  points c+s*dc and c+0.5*s*dc, with dc = the proposed change in the vector
  of parameters, and with step length s initially at 1. s is halved until
  minus the objective function stops declining.
@
proc (2) = squeeze(s_c,s_dc);
local s_c1,s_lm,s_itr,lc1,s_c2,lc2,s_f1,s_f2;


    s_c1=s_c - s_dc; s_lm=1/2; s_itr=1;
    s_f1 = deff(s_c1);
    lc1 = s_f1'w*s_f1;
    clear s_f1;

  do until s_itr > maxsqez;

    s_c2=s_c-s_lm*s_dc;
    s_f2 = deff(s_c2);
    lc2 = s_f2'w*s_f2;
    clear s_f2;

    if lc1 <= lc2 and lc1 <= obj;

       retp(s_c1,s_itr-1); goto eoproc;

    else;

       s_c1=s_c2; s_lm=s_lm/2; s_c2=s_c - s_lm*s_dc; lc1=lc2;
       s_itr=s_itr+1;

    endif;

  endo;

retp(s_c2,s_itr-1);
eoproc:

endp;

@  ------------------------------------------------------ @
@ Subroutine to compute optimal weighting matrix. @

proc optw(a);

local f,ff;

f = zeros(n,neq);

f[.,1] = y2x - a[1,1] + (-2*Eyx_z'inEzz*zy_d' - Ey2_z'inEzz*zx_d')';
f[.,2] = yx2 - a[2,1] + (-Ex2_z'inEzz*zy_d' - 2*Eyx_z'inEzz*zx_d')';

ff = moment(f,0);

clear f;

    retp(invpd(ff./n));
    clear ff;

endp;

eop:
end;   @ This ends the program and closes all files.   @
@ ----------------------- END OF PROGRAM ------------------------------- @
