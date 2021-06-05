new ,60000;
ttt = time;
ddd = date;

@------Indicate whether the last two observations are to be used. -@

   outlier = 0;

@------Name of the Gauss data set.---------------------------------@

dataset1 = "ewdset";

outwidth  132;


@-------Pick the q proxy.------------------------------------------@

qqq = 1;

qidx = 5;
output file = muxmuy.out reset;

@--------Open GAUSS data set. -------------------------------------@

    closeall f1;
    open f1=^dataset1;

@------------------Set the number of years.-----------------------@

    nyr = 4;


@-------These loops do all possible combinations of models---------@

@-----------------Set the number of perfectly measured regressors---@

    nz = 1;

@-----------------Set the maximum number of iterations.------------@

maxiter = 299;

@-----------------Set the maximum number of rounds.------------@

rlim = 1;

@-----------------Set the limit on the number of squeezes.---------@

    maxsqez=240;

@-------Set the flag for the standard error calculations.----------@

bleh = 2;

/* bleh = 1 acts as if Ey and Ex known.  bleh = 2 optimally accounts for the
   fact that they are not known. */

@-----------------Set the number of observations.------------------@

    nco = rowsf(f1)/nyr;


cls;



   year = 1;
   do until year > nyr;


    @ --------- Read in the data and define the variables. ----------@

      dta=readr(f1,nco);


    n=rows(dta);


@ cusip sic year ik q cf k da rate payrat sales btm @
@ 1     2   3    4  5 6  7 8  9    10     11    12  @


@------------This part has to be changed by hand to alter the number of regressors.----@



    ww2 = dta[.,9] .== 9999;

    for ii(0,1,1);


    y = selif(dta[.,4],ww2-ii);
    z = selif(dta[.,6],ww2-ii);
    x = selif(dta[.,5],ww2-ii);



    iz = ones(rows(z),1)~z;


    mux = (invpd(moment(iz,0)))* (iz'x);
    muy = (invpd(moment(iz,0)))* (iz'y);
    varz = meanc((z-meanc(z))^2);


    if ii==1;
    year+1991;; "-------------------Unconstrained--------------";;   "muy    ";; muy[2,1];
    year+1991;; "-------------------Unconstrained--------------";;   "mux    ";; mux[2,1];
    year+1991;; "-------------------Unconstrained--------------";;   "ratio  ";; muy[2,1]/mux[2,1];
    year+1991;; "-------------------Unconstrained--------------";;   "varz   ";; varz;
    else;
    year+1991;; "-------------------Constrained----------------";;   "muy    ";; muy[2,1];
    year+1991;; "-------------------Constrained----------------";;   "mux    ";; mux[2,1];
    year+1991;; "-------------------Constrained----------------";;   "ratio  ";; muy[2,1]/mux[2,1];
    year+1991;; "-------------------Constrained----------------";;   "varz   ";; varz;
    endif;









    endfor;

    year = year + 1;

endo;     @This ends the trial loop.@

"Start time:  ";; timestr(ttt);
"End time:    ";; timestr(0);
"Start date:  ";; datestr(ddd);
"End date:    ";; datestr(0);

output off;
end;
