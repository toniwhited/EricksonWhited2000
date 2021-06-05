new ,60000;
ttt = time;
ddd = date;

@------Indicate whether the last two observations are to be used. -@

   outlier = 0;

@------Name of the Gauss data set.---------------------------------@

dataset1 = "ewdset";

outwidth  132;

qidx = 5;
output file = ewgmm.out reset;

@--------Open GAUSS data set. -------------------------------------@

    closeall f1;
    open f1=^dataset1;

@------------------Set the number of years.-----------------------@

    nyr = 4;


@-------These loops do all possible combinations of models---------@

   strup = 1;
   do while strup <=1;

   ndum = 2;
   do while ndum >= 1;

   pick = 1;
   do while pick <= (3!)/((ndum!)*((3-ndum)!)) + 1*(ndum==1);

     call seekr(f1,1);

     if strup == 1; "With Straight-Up";
     else;          "Without Straight-Up";
     endif;
     if ndum == 3;
                       "All Three Dummies";            pickidx = { 1 2 3 };
     elseif ndum == 2;
           if pick==1; "Dividend and Bond Dummies";    pickidx = { 1 2 };
       elseif pick==2; "Dividend and Size Dummies";    pickidx = { 1 3 };
       elseif pick==3; "Bond and Size Dummies";        pickidx = { 2 3 };
       endif;
     elseif ndum == 1;
           if pick==1; "Dividend Dummy";               pickidx = { 1 };
       elseif pick==2; "Bond Dummy";                   pickidx = { 2 };
       elseif pick==3; "Size Dummy";                   pickidx = { 3 };
       elseif pick==4; "Both Dummy";                   pickidx = { 2 3 };

       endif;
     endif;


@-----------------Set the number of perfectly measured regressors---@

    nz = 1+(strup+1)*ndum;

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

@-----------This starts the year loop.--------------------------------@

   nestim = 3;

@-------------------Parameters-----------------------@


@OLS@

   isave   = zeros(2,nyr);         @ Intercept                                     @
   bsave   = zeros(2,nyr);         @ Coefficient on chi                            @
   zsave   = zeros(2*nz,nyr);      @ Coefficients on perfectly measured regressors @
   zisave  = zeros(2,nyr);         @ Sum of the interaction terms                  @
   zdsave  = zeros(2,nyr);         @ Sum of the straight up dummies                @
   rbsave  = zeros(2,nyr);         @ Reverse chi                                   @
   rzsave  = zeros(2*nz,nyr);      @ Reverse perfectly measured regressors         @
   rsave   = zeros(2,nyr);         @ R2                                            @

@GMM@

   asave   = zeros(2,nyr*nestim);    @ Intercept                                     @
   csave   = zeros(2,nyr*nestim);    @ Coefficient on chi                            @
   dsave   = zeros(2*nz,nyr*nestim); @ Coefficients on perfectly measured regressors @
   disave  = zeros(2,nyr*nestim);    @ Sum of the interaction terms                  @
   ddsave  = zeros(2,nyr*nestim);    @ Sum of the straight up dummies                @
   tsave   = zeros(2,nyr*nestim);    @ tau2                                          @
   psave   = zeros(2,nyr*nestim);    @ rho2                                          @
   chisave = zeros(2,nyr*nestim);    @ J-statistic                                   @

@-------------------Minimum Distance Estimates-------@

@GMM@

   casave  = zeros(4,nestim);     @ Intercept                                     @
   ccsave  = zeros(4,nestim);     @ Coefficient on chi                            @
   cdsave  = zeros(4*nz,nestim);  @ Coefficients on perfectly measured regressors @
   ctsave  = zeros(4,nestim);     @ tau2                                          @
   cpsave  = zeros(4,nestim);     @ rho2                                          @
   cdisave = zeros(4,nestim);     @ Sum of the interaction terms                  @
   cddsave = zeros(4,nestim);     @ Sum of the straight up dummies                @

@OLS@

   bisave  = zeros(4,1);          @ Intercept                                     @
   bbsave  = zeros(4,1);          @ Coefficient on chi                            @
   brsave  = zeros(4,1);          @ R2                                            @
   bzsave  = zeros(4*nz,1);       @ Coefficients on perfectly measured regressors @
   bzisave = zeros(4,1);          @ Sum of the interaction terms                  @
   bzdsave = zeros(4,1);          @ Sum of the straight up dummies                @



@-------------------Influence Functions--------------@

@GMM@

   iasave  = zeros(nco,nyr*nestim);
   icsave  = zeros(nco,nyr*nestim);
   idsave  = zeros(nco*nz,nyr*nestim);
   ipsave  = zeros(nco,nyr*nestim);
   itsave  = zeros(nco,nyr*nestim);

@OLS@

   iisave  = zeros(nco,nyr);
   ibsave  = zeros(nco,nyr);
   irsave  = zeros(nco,nyr);
   izsave  = zeros(nco*nz,nyr);

cls;


call seekr(f1,1);

   year = 1;
   do until year > nyr;


    @ --------- Read in the data and define the variables. ----------@

      dta=readr(f1,nco);


    n=rows(dta);


@ cusip sic year ik q cf k da rate payrat sales btm @
@ 1     2   3    4  5 6  7 8  9    10     11    12  @


@------------This part has to be changed by hand to alter the number of regressors.----@


    y = dta[.,4];
    z = dta[.,6];
    x = dta[.,qidx];

    ww = dta[.,10];
    ww = ww .== 0;
    ww1 = 1 - ww;

    brate = sortc(dta[.,9],1);
    ww2 = dta[.,9] .== 9999;

    wz1 = ww1.*z;
    wz2 = ww2.*z;



both=1;
if both==0;
  capital=0;
  if capital==0;

    big = sortc(dta[.,15],1);
    one3b = big[floor(1*rows(big)/3),1];
    two3b = big[floor(2*rows(big)/3),1];

    ww3 = dta[.,15] .< one3b;
    ww4 = dta[.,15] .> two3b;
  elseif capital==1;

    big = sortc(dta[.,7],1);
    one3b = big[floor(1*rows(big)/3),1];
    two3b = big[floor(2*rows(big)/3),1];

    ww3 = dta[.,7] .< one3b;
    ww4 = dta[.,7] .> two3b;
  endif;


else;

  big1   = sortc(dta[.,7],1);
  one3b1 = big1[floor(1*rows(big1)/3),1];
  two3b1 = big1[floor(2*rows(big1)/3),1];
  big2   = sortc(dta[.,15],1);
  one3b2 = big2[floor(1*rows(big2)/3),1];
  two3b2 = big2[floor(2*rows(big2)/3),1];

  ww3 = dta[.,7] .< one3b1 .and dta[.,15] .< one3b2;
  ww4 = dta[.,7] .> two3b1 .and dta[.,15] .> two3b2;


endif;
    wz3 = ww3.*z;
    wz4 = ww4.*z;


    allwz = wz1~wz2~wz3;
    allww = ww1~ww2~ww3;

    thiswz = allwz[.,pickidx];
    thisww = allww[.,pickidx];

    icept=ones(n,1);
      if strup == 0;
        za = z~thiswz;
      elseif strup == 1;
        za = z~thiswz~thisww;
        if ndum==1 and pick==4;
        addon = prodc(allww[.,pickidx]');
        zaddon = addon.*z;
        za = z~zaddon~addon;
        endif;
      endif;

    iz = icept~za;

    mux = (invpd(moment(iz,0)))* (iz'x);
    muy = (invpd(moment(iz,0)))* (iz'y);
    y_d = y - iz*muy;
    x_d = x - iz*mux;

    y2_d = y_d^2;
    yx_d = (x_d .* y_d);
    x2_d = x_d^2;
    y2x_d = (x_d .* y2_d);
    yx2_d = (x2_d .* y_d);
    y3x_d = y2x_d .* y_d;
    y2x2_d = y2_d .* x2_d;
    yx3_d = yx2_d .* x_d;

    x3_d = x_d.*x2_d;
    y3_d = y_d.*y2_d;
    y4x_d = y3_d.*yx_d;
    y3x2_d = y3_d.*x2_d;
    y2x3_d = y2_d.*x3_d;
    yx4_d = yx_d.*x3_d;
    y4_d = y2_d^2;
    x4_d = x2_d^2;
    y5x_d = y_d.*y4x_d;
    y4x2_d = y_d.*y3x2_d;
    y3x3_d = y_d.*y2x3_d;
    y2x4_d = y_d.*yx4_d;
    yx5_d = yx4_d.*x_d;

    y5_d = y2_d.*y3_d;
    x5_d = x2_d.*x3_d;
    y6x_d = y5_d.*yx_d;
    y5x2_d = y5_d.*x2_d;
    y4x3_d = y4_d.*x3_d;
    y3x4_d = y3_d.*x4_d;
    y2x5_d = y2_d.*x5_d;
    yx6_d = yx_d.*x5_d;

/* VERY IMPORTANT:  The following variables MUST be in the same order that
                    they appear in the subroutine deff.
*/

    mom = { y2_d yx_d x2_d y2x_d yx2_d y3x_d y2x2_d yx3_d x3_d y3_d y4x_d
            y3x2_d y2x3_d yx4_d y4_d x4_d y5x_d y4x2_d y3x3_d y2x4_d yx5_d
            y5_d x5_d y6x_d y5x2_d y4x3_d y3x4_d y2x5_d yx6_d };

    Ey2_d = meanc(y2_d);
    Eyx_d = meanc(yx_d);
    Ex2_d = meanc(x2_d);
    Ey2x_d = meanc(y2x_d);
    Eyx2_d = meanc(yx2_d);
    Ey3x_d = meanc(y3x_d);
    Ey2x2_d = meanc(y2x2_d);
    Eyx3_d = meanc(yx3_d);

    Ex3_d     = meanc(x3_d);
    Ey3_d     = meanc(y3_d);
    Ey4x_d    = meanc(y4x_d);
    Ey3x2_d   = meanc(y3x2_d);
    Ey2x3_d   = meanc(y2x3_d);
    Eyx4_d    = meanc(yx4_d);
    Ey4_d     = meanc(y4_d);
    Ex4_d     = meanc(x4_d);
    Ey5x_d    = meanc(y5x_d);
    Ey4x2_d   = meanc(y4x2_d);
    Ey3x3_d   = meanc(y3x3_d);
    Ey2x4_d   = meanc(y2x4_d);
    Eyx5_d    = meanc(yx5_d);
    Ey5_d     = meanc(y5_d);
    Ex5_d     = meanc(x5_d);
    Ey6x_d    = meanc(y6x_d);
    Ey5x2_d   = meanc(y5x2_d);
    Ey4x3_d   = meanc(y4x3_d);
    Ey3x4_d   = meanc(y3x4_d);
    Ey2x5_d   = meanc(y2x5_d);
    Eyx6_d    = meanc(yx6_d);

    Emom = { Ey2_d Eyx_d Ex2_d Ey2x_d Eyx2_d Ey3x_d Ey2x2_d Eyx3_d Ex3_d
             Ey3_d Ey4x_d Ey3x2_d Ey2x3_d Eyx4_d Ey4_d Ex4_d Ey5x_d Ey4x2_d
             Ey3x3_d Ey2x4_d Eyx5_d Ey5_d Ex5_d Ey6x_d Ey5x2_d Ey4x3_d
             Ey3x4_d Ey2x5_d Eyx6_d };

/*  Make the moments for the second stage of the partialling estimator. */


    Ezx  = iz'x/n;

/* Make the moments for the rest of the summary statistics. */

   Ey  = meanc(y);
   Ex  = meanc(x);
   Ez  = meanc(z);
   Ey2 = meanc(y^2);
   Ex2 = meanc(x^2);
   Ez2 = meanc(z^2);


/* Make the moments for the correct standard errors for beta. */



    y6_d  = y_d^6;
    x6_d  = x_d^6;

    Ey2_z   =   meanc(y2_d   .* iz);
    Eyx_z   =   meanc(yx_d   .* iz);
    Ex2_z   =   meanc(x2_d   .* iz);
    Ey2x_z  =   meanc(y2x_d  .* iz);
    Eyx2_z  =   meanc(yx2_d  .* iz);
    Ey3x_z  =   meanc(y3x_d  .* iz);
    Ey2x2_z =   meanc(y2x2_d .* iz);
    Eyx3_z  =   meanc(yx3_d  .* iz);
    Ey3_z   =   meanc(y3_d   .* iz);
    Ex3_z   =   meanc(x3_d   .* iz);
    Ey4x_z    = meanc(y4x_d.*iz);
    Ey3x2_z   = meanc(y3x2_d.*iz);
    Ey2x3_z   = meanc(y2x3_d.*iz);
    Eyx4_z    = meanc(yx4_d.*iz);
    Ey4_z     = meanc(y4_d.*iz);
    Ex4_z     = meanc(x4_d.*iz);
    Ey5x_z    = meanc(y5x_d.*iz);
    Ey4x2_z   = meanc(y4x2_d.*iz);
    Ey3x3_z   = meanc(y3x3_d.*iz);
    Ey2x4_z   = meanc(y2x4_d.*iz);
    Eyx5_z    = meanc(yx5_d.*iz);
    Ey5_z     = meanc(y5_d.*iz);
    Ex5_z     = meanc(x5_d.*iz);
    Ey6x_z    = meanc(y6x_d.*iz);
    Ey5x2_z   = meanc(y5x2_d.*iz);
    Ey4x3_z   = meanc(y4x3_d.*iz);
    Ey3x4_z   = meanc(y3x4_d.*iz);
    Ey2x5_z   = meanc(y2x5_d.*iz);
    Eyx6_z    = meanc(yx6_d.*iz);
    Ey6_z     = meanc(y6_d.*iz);
    Ex6_z     = meanc(x6_d.*iz);

    Ezz     =   moment(iz,0)/n;
    inEzz   =   invpd(Ezz);
    zy_d    =   iz.*y_d;
    zx_d    =   iz.*x_d;

    Eza   = meanc(za)';

/* The OLS estimator: */

    des = icept~x~za;
    dd = invpd(moment(des,0));
    allb = dd*des'y;
    isave[1,year] = allb[1,1];
    bsave[1,year] = allb[2,1];
    for qq(1, nz, 1); zsave[(2*qq-1),year] = allb[2+qq,1]; endfor;
    uhat = y - des*allb;
    wh = uhat.*des;
    seb = sqrt(diag(dd*(wh'wh)*dd));
    naiveseb = sqrt(diag(dd*meanc(uhat^2)));
    isave[2,year] = seb[1,1];
    bsave[2,year] = seb[2,1];
    for qq(1, nz, 1); zsave[(2*qq),year] = seb[2+qq,1]; endfor;

    rsave[1,year] = 1 - moment(uhat,0)/(moment((y - Ey),0));

    olsinflnc = invpd(moment(des,0)/n)*((des.*uhat)');
    iisave[.,year] = olsinflnc[1,.]';
    ibsave[.,year] = olsinflnc[2,.]';
    for qq(1, nz, 1); izsave[((qq-1)*n+1):(qq*n),year] = olsinflnc[2+qq,.]'; endfor;


@----------First make the influence function for sigma_xz.--------------------@
if nz > 0;

xza = x~za;
Exza = meanc(xza);
nreg=nz+1;
sigxz = moment((xza)-Exza',0)/n;


vecsigxz = zeros(nreg+nreg*(nreg-1)/2,1);
vecsigxz[1:nreg,1] = diag(sigxz);

counter =  nreg+nreg*(nreg-1)/2;
for qq(nreg, 2, -1);
  for ee(qq-1, 1, -1);
    vecsigxz[counter,1] = sigxz[ee,qq];
    counter=counter-1;
  endfor;
endfor;


phixz = zeros(n,nreg+nreg*(nreg-1)/2);
phixz[.,1:nreg] = ((xza) - Exza').*((xza) - Exza');

counter =  nreg+nreg*(nreg-1)/2;
for qq(nreg, 2, -1);
  for ee(qq-1, 1, -1);
    phixz[.,counter] = (xza[.,ee] - Exza[ee,1]).*(xza[.,qq] - Exza[qq,1]);
    counter=counter-1;
  endfor;
endfor;

sigy = moment((y - Ey),0)/n;
phiy = (y - Ey)^2 - sigy;

bigphi = (olsinflnc[2:rows(olsinflnc),.])|(phixz')|(phiy');

@--------------Now make the derivative matrix.-------------------------------@

gee = zeros(nreg+(nreg+nreg*(nreg-1)/2)+1,1);


@-------------First, the derivative wrt the OLS coefficients.----------@

for qq(1, nreg, 1);

gee[qq,1] = (2*allb[2:nreg+1,1]'sigxz[.,qq])/sigy;

endfor;

@------------derivatives wrt the first part of sigxz--@

for qq(1, nreg, 1);

gee[nreg+qq,1] = (allb[qq+1,1]^2)/sigy;

endfor;

@------------derivatives wrt the second part of sigxz--@

counter = nreg+(nreg+nreg*(nreg-1)/2);
for qq(nreg, 2, -1);
  for ee(qq-1, 1, -1);
  gee[counter,1] = 2*allb[ee+1,1]*allb[qq+1,1]/sigy;
  counter=counter-1;
  endfor;
endfor;

@------------derivatives wrt sigy--@

counter = nreg+(nreg+nreg*(nreg-1)/2);

gee[counter+1,1] = -rsave[1,year]/sigy;

@-----------Save the influence function.-------------@

else;

@----------First make the influence function for sigma_xz.--------------------@
sigx = Ex2-Ex^2;
phix = (x - Ex)^2 - sigx;

sigy = moment((y - Ey),0)/n;
phiy = (y - Ey)^2 - sigy;

bigphi = (olsinflnc[2:rows(olsinflnc),.])|(phix')|(phiy');

@--------------Now make the derivative matrix.--------------------------------@

gee = zeros(3,1);

@-------------First, the derivative wrt the OLS coefficients.----------@

gee[1,1] = (2*allb[2,1]*sigx)/sigy;

@------------derivatives wrt the first part of sigxz--@

gee[2,1] = (allb[2,1]^2)/sigy;

@------------derivatives wrt sigy--@

gee[3,1] = -rsave[1,year]/sigy;

@-----------Save the influence function.-------------@

endif;

irsave[1:n,year] = -bigphi'gee;
 rsave[2,year] =  sqrt(moment(irsave[.,year],0)/n^2);

@-------------------This part does the test that alpha_1 plus alpha_2 = 0----------@

    gee  = des'uhat/n;
    capg = -(des'des)/n;
    omega = wh'wh/n;

    @---------------Test that the noninteracted dummies are jointly significant.--@
if strup == 1;
       capA = zeros(ndum,(3+ndum))~eye(ndum);
       wald = n*allb'capA'invpd(capA*invpd(capG'invpd(omega)*capG)*capA')*capA*allb;
@      "OLS Test for joint significance of the dummies";                        @
@      "Wald Statistic " wald;                                                  @
@      "P-value        " cdfchic(wald,3);                                       @
endif;

    des = icept~y~za;
    dd = invpd(moment(des,0));
    allbr = dd*des'x;
    rbsave[1,year]  = 1/allbr[2,1];
    for qq(1, nz, 1); rzsave[(2*qq-1),year] = -allbr[qq+2,1]/allbr[2,1]; endfor;
    uhat = x - des*allbr;
    wh = uhat.*des;
    vc=(dd*(wh'wh)*dd);
     gee = zeros(nz+2,nz+2);

     gee[1,1] = -1/allbr[2,1];
     gee[1,2] = allbr[1,1]/(allbr[2,1]^2);
     gee[2,2] = -allbr[2,1]^(-2);
     for qq(1, nz, 1); gee[qq+2,qq+2] = -1/allbr[2,1]; endfor;
     for qq(1, nz, 1); gee[qq+2,2] = allbr[qq+2,1]/allbr[2,1]; endfor;

     vcr2t = gee*vc*gee';
     serb = sqrt(diag(vcr2t));
     rbsave[2,year] = serb[2,1];
     for qq(1, nz, 1); rzsave[(2*qq),year] = serb[qq+1,1]; endfor;

/* The GMM estimators: */

/* bleh = 1 acts as if Ey and Ex known.  bleh = 2 optimally accounts for the
   fact that they are not known.
*/
    estim = 1;
    do while estim <= nestim;

@--Once exaclly identified and four iterative estimators. -----------------@

if estim==1;
 nc = 5;
 neq = 5;
elseif estim==2;
 nc = 6;
 neq = 8;
elseif estim==3;
 nc = 9;
 neq = 14;
elseif estim==4;
 nc = 12;
 neq = 21;
elseif estim==5;
 nc = 15;
 neq = 29;
endif;
    @---------- Here we input the starting values.----------------@
    c = zeros(nc,1);

    c[1,1] = Ey2x_d/Eyx2_d;

@   if estim>1;c[1,1]=bsave[1,year];endif;@


    c[2,1] = Eyx_d/c[1,1];
    c[3,1] = Ey2_d - (c[1,1]^2)*c[2,1];
    c[4,1] = Ex2_d - c[2,1];
    c[5,1] = Eyx2_d/c[1,1];
if estim > 1;
    c[6,1] = (Eyx3_d - 3*c[1,1]*c[2,1]*c[4,1])/c[1,1];
endif;
if estim > 2;
    c[7,1] = Ex3_d - c[5,1];
    c[8,1] = Ey3_d - (c[1,1]^3)*c[5,1];
    c[9,1] = (Eyx4_d/c[1,1]) - 6*c[5,1]*c[4,1] - 4*c[2,1]*c[7,1];
endif;
if estim > 3;
    c[10,1] = Ey4_d - (c[1,1]^4)*c[6,1] - 6*(c[1,1]^2)*c[2,1]*c[3,1];
    c[11,1] = Ex4_d - c[6,1] - 6*c[2,1]*c[4,1];
    c[12,1] = (Eyx5_d/c[1,1]) - 10*c[6,1]*c[4,1] - 10*c[5,1]*c[7,1]
              - 5*c[2,1]*c[11,1];
endif;
if estim > 4;
    c[13,1] = Ey5_d - (c[1,1]^5)*c[9,1] - 10*(c[1,1]^3)*c[5,1]*c[3,1]
              - 10*(c[1,1]^2)*c[2,1]*c[8,1];
    c[14,1] = Ex5_d - c[9,1] - 10*c[5,1]*c[4,1] - 10*c[2,1]*c[7,1];
    c[15,1] = (Eyx6_d/c[1,1]) - 15*c[9,1]*c[4,1] - 20*c[6,1]*c[7,1]
              - 15*c[5,1]*c[11,1] - 6*c[2,1]*c[14,1];
endif;

if estim>1;

     @ --------The program creates the weighting matrix. ------------@

  ff = optw(c,mom,Emom,bleh);
  w = invit(moment(ff,0)./n);

     @ ------------- Start of iteration loop ------------------------@
iter = 1;
    tol=1e-9;                 @ Set the convergence criterion. @
    dc=1;                     @ Initialize the step length. @
do until abs(dc) < tol;

    f=deff(c,Emom,0);         @ The program jumps to the subroutine that
                               defines the f vector. @

    if estim == 2;
    g=grad2(c);                @ The program jumps to the subroutine that
                                computes analytic derivatives. The matrix
                                of partials of f with respect to the
                                parameters is called g. @
    elseif estim == 3;
    g=grad3(c);
    elseif estim == 4;
    g=grad4(c);
    elseif estim == 5;
    g=grad5(c);
    endif;

    obj = f'w*f;              @ This computes the value of the objective
                                function.@

    gwg= g'w*g;               @ This uses the GAUSS-NEWTON method
                                to compute full step dc. @
    gwf = g'w*f;

    dc = solvit(gwf,gwg);

 if maxsqez > 0;              @ This jumps to the subroutine that
                                adjusts the step length. @
   { c_new,sqz } = squeeze(c,dc);
 else;
   c_new=c - dc;
 endif;


 dc=c_new-c;                     @ Update variables for the next iteration. @
 c=c_new;
/*
/*--------------Print out the results for the current iteration------------*/
cls;
"Number of Squeezes = ";; sqz;
"Objective Function = ";; obj*n;

"            Value         Step";
 format 14,6;
 mm=1;
 do until mm > nc;
  $"  ";; c[mm,1];; "  ";; dc[mm,1];
  mm=mm+1;
  endo;
/*-----------------End of print out of current iteration--------------*/
*/
 iter=iter +1;
 if iter >= maxiter;                     @ Quit iterating if necessary. @
   goto escape;
 endif;
endo;
@ ------ End of iteration loop ------------------------------------ @
escape:
endif;

@ Compute t-ratios, etc.                                            @

  ff = optw(c,mom,Emom,bleh);
  f = deff(c,Emom,0);

  w = invit(moment(ff,0)./n);

    if estim == 1;
    g=grad1(c);
    elseif estim == 2;
    g=grad2(c);
    elseif estim == 3;
    g=grad3(c);
    elseif estim == 4;
    g=grad4(c);
    elseif estim == 5;
    g=grad5(c);
    endif;

  gwg = g'w*g;
  vc=invit(gwg)/n;
    stderr=sqrt(diag(vc));
    inflnc=-n*vc*g'w*ff';

    csave[1,nestim*(year-1)+estim] = c[1,1];
    csave[2,nestim*(year-1)+estim] = stderr[1,1];
    icsave[.,nestim*(year-1)+estim] = inflnc[1,.]';

if estim > 1;
chisq = obj*n;
chisave[1,nestim*(year-1)+estim] = chisq;
degf=neq - nc;
chisave[2,nestim*(year-1)+estim] = cdfchic(chisq,degf);
endif;

/* Make Delta */

c11 = c[1,1];
for qq(1, nz, 1); dsave[(qq*2-1),nestim*(year-1)+estim] = muy[qq+1,1] - c11*mux[qq+1,1]; endfor;


gee = zeros(((nz+1)*2+1),nz);

gee[2:(nz+1),1:nz] = eye(nz);
gee[(nz+3):((nz+1)*2),1:nz] = -c11*eye(nz);

for qq(1, nz, 1); gee[((nz+1)*2+1),qq] = -mux[qq+1,1]; endfor;

@---Correct standard error and influence function---@

bigphi = ( (inEzz*zy_d')|(inEzz*zx_d') )|(-inflnc[1,.]);
avar = moment(bigphi',0)./n^2;
dstd = sqrt(diag(gee'avar*gee));
omega = gee'avar*gee;
phidel = -bigphi'gee;

for qq(1, nz, 1); dsave[qq*2,nestim*(year-1)+estim] = dstd[qq,1]; endfor;

for qq(1, nz, 1); idsave[((qq-1)*n+1):(qq*n),nestim*(year-1)+estim] = phidel[.,qq]; endfor;


/*  Make the stupid intercept */

asave[1,nestim*(year-1)+estim] = muy[1,1] - c11*mux[1,1];

gee = zeros(((nz+1)*2+1),1);

gee[1,1] = 1;
gee[(nz+2),1] = -c11;
gee[((nz+1)*2+1),1] = -mux[1,1];

@---Correct standard error and influence function---@

bigphi = ( (inEzz*zy_d')|(inEzz*zx_d') )|(-inflnc[1,.]);
avar = moment(bigphi',0)./n^2;
aastd = sqrt(diag(gee'avar*gee));
phidela = -bigphi'gee;

asave[2,nestim*(year-1)+estim] = aastd[1,1];
iasave[1:n,nestim*(year-1)+estim] = phidela[.,1];


@----------First make the influence function for sigma_z.----------------------@

sigz = moment(za-meanc(za)',0)/n;


vecsigz = zeros(nz+nz*(nz-1)/2,1);
vecsigz[1:nz,1] = diag(sigz);

counter =  nz+nz*(nz-1)/2;
for qq(nz, 2, -1);
  for ee(qq-1, 1, -1);
    vecsigz[counter,1] = sigz[ee,qq];
    counter=counter-1;
  endfor;
endfor;


phiz = zeros(n,nz+nz*(nz-1)/2);
phiz[.,1:nz] = (za - Eza).*(za - Eza);

counter =  nz+nz*(nz-1)/2;
for qq(nz, 2, -1);
  for ee(qq-1, 1, -1);
    vecsigz[counter,1] = sigz[ee,qq];
    phiz[.,counter] = (za[.,ee] - Eza[1,ee]).*(za[.,qq] - Eza[1,qq]);
    counter=counter-1;
  endfor;
endfor;


phiz = phiz - vecsigz';

phimuy = inEzz*zy_d';
phimux = inEzz*zx_d';

numer = (muy[2:nz+1,1]'sigz*muy[2:nz+1,1] + c[1,1]^2*c[2,1]);
rho = numer/(numer + c[3,1]);

numer = (mux[2:nz+1,1]'sigz*mux[2:nz+1,1] + c[2,1]);
tau = numer/(numer + c[4,1]);

@---------Make the influence functions for the standard errors for rho2 and tau2.@

bigphi = phimux[2:nz+1,.]|phimuy[2:nz+1,.]|(phiz')|(-inflnc[1:4,.]);
avar = moment(bigphi',0)./n^2;

gee = zeros(2*nz+(nz+nz*(nz-1)/2)+4,2); @ First column is for rho2 and the second is for tau2. @


@-----------First, do rho2---------------------------@


numer = muy[2:nz+1,1]'sigz*muy[2:nz+1,1]+c[1,1]^2*c[2,1];
denom =  numer+c[3,1]; ;

@------------derivatives wrt muy---------------------@

for qq(1, nz, 1);

gee[nz+qq,1] = (2*muy[2:nz+1,1]'sigz[.,qq])/denom - numer*(2*muy[2:nz+1,1]'sigz[.,qq])/(denom^2);

endfor;
@------------derivatives wrt the first part of sigz--@

for qq(1, nz, 1);

gee[2*nz+qq,1] = (muy[qq+1,1]^2)/denom - numer/(denom^2)*muy[qq+1,1]^2;

endfor;

@------------derivatives wrt the second part of sigz--@

counter = 2*nz+(nz+nz*(nz-1)/2);
for qq(nz, 2, -1);
  for ee(qq-1, 1, -1);
  gee[counter,1] = 2*muy[ee+1,1]*muy[qq+1,1]/denom - 2*numer/(denom^2)*muy[ee+1,1]*muy[qq+1,1];
  counter=counter-1;
  endfor;
endfor;

@------------derivatives wrt  c--@

counter = 2*nz+(nz+nz*(nz-1)/2);

gee[counter+1,1] = 2*c[1,1]*c[2,1]/denom - 2*numer/(denom^2)*c[1,1]*c[2,1];

gee[counter+2,1] = (c[1,1]^2)/denom - numer/(denom^2)*c[1,1]^2;

gee[counter+3,1] = -numer/(denom^2);


@-----------Now for Tau2------------------------@


numer = mux[2:nz+1,1]'sigz*mux[2:nz+1,1]+c[2,1];
denom = numer+c[4,1]; ;



@------------derivatives wrt mux---------------------@

for qq(1, nz, 1);

gee[qq,2] = (2*mux[2:nz+1]'sigz[.,qq])/denom - numer*(2*mux[2:nz+1]'sigz[.,qq])/(denom^2);

endfor;

@------------derivatives wrt the first part of sigz--@

for qq(1, nz, 1);

gee[2*nz+qq,2] = (mux[qq+1,1]^2)/denom - numer/(denom^2)*mux[qq+1,1]^2;

endfor;

@------------derivatives wrt the second part of sigz--@

counter = 2*nz+(nz+nz*(nz-1)/2);
for qq(nz, 2, -1);
  for ee(qq-1, 1, -1);
  gee[counter,2] = 2*mux[ee+1,1]*mux[qq+1,1]/denom - 2*numer/(denom^2)*mux[ee+1,1]*mux[qq+1,1];
  counter=counter-1;
  endfor;
endfor;

gee[2*nz+(nz+nz*(nz-1)/2)+4,2] = -numer/(denom^2);


vcrhotau = gee'avar*gee;

phirho = -bigphi'gee[.,1];
phitau = -bigphi'gee[.,2];

ipsave[.,nestim*(year-1)+estim] = phirho;
itsave[.,nestim*(year-1)+estim] = phitau;


/*  Now save both things. */

    psave[1,nestim*(year-1)+estim]    = rho;
    tsave[1,nestim*(year-1)+estim]    = tau;

    psave[2,nestim*(year-1)+estim]    = sqrt(vcrhotau[1,1]);
    tsave[2,nestim*(year-1)+estim]    = sqrt(vcrhotau[2,2]);

/*
@     "Year = " year+1991;                                                              @
@     "Estimator = " estim+2;                                                           @

    @-----------Test that beta + alpha_3 = 0----------------------@

     capA = { 1 1 };
     sigmahat = sqrt(capA*omega*capA'/n);

@   "Wald test for (alpha_1+alpha_2=0)        " ((idel+ddel)/sigmahat)^2/n;             @
@   "P-value                                  " cdfchic(((idel+ddel)/sigmahat)^2/n,1);  @
@   "Standard Error for (alpha_1+alpha_2)     " sigmahat;                               @
*/


    @-----------Test that the noninteractive dummies are jointly significant---------------@
if strup == 1;
     capA = zeros(ndum,ndum+1)~eye(ndum);
     idx = seqa(1,2,rows(dsave)/2);
     wald =
     dsave[idx,nestim*(year-1)+estim]'capA'invpd(capA*(omega)*capA')*capA*dsave[idx,nestim*(year-1)+estim];
endif;

@    "GMM Test for joint significance of the dummies";                                  @
@    "Estimator = " estim+2;                                                            @
@    "Wald Statistic " wald;                                                            @
@    "P-value        " cdfchic(wald,3);                                                 @


    estim = estim + 1;
    endo;                       @ The end of the estimator loop. @

    year = year + 1;

endo;     @This ends the trial loop.@



/* This part does the classical minimum distance estimation of beta and alpha_1. */

param = 1; @ 1 is beta, 2 is rho2, 3 is tau2,
             4 to nz+3 are the alphas, and the end is the intercept. @
do while param <= nz+4;

  idx = seqa(1,nestim,nyr) - 1;

estim = 1;
do while estim <= nestim;

  idx = idx+1;

  if param == 1;
    saveit = csave;
    w = invit(moment(icsave[.,idx],0)/n);
  elseif param == 2;
    saveit = psave;
    w = invit(moment(ipsave[.,idx],0)/n);
  elseif param == 3;
    saveit = tsave;
    w = invit(moment(itsave[.,idx],0)/n);
  elseif param == nz+4;
    saveit = asave;
    w = invit(moment(iasave[.,idx],0)/n);
  endif;

  for qq(1, nz, 1);
  if param == qq+3;
    saveit = dsave[qq*2-1,.];
    w = invit(moment(idsave[((qq-1)*n+1):(qq*n),idx],0)/n);
  endif;
  endfor;


    g=ones(nyr,1);
    theta = invpd(g'w*g)*g'w*saveit[1,idx]';

    f=saveit[1,idx]'-theta;
    stderror = sqrt(invpd(g'w*g)/n);
    cmdtest = n*f'w*f;

    if param == 1;      ccsave[1,estim] = theta; ccsave[2,estim] = cmdtest;  ccsave[3,estim] = cdfchic(cmdtest,3); ccsave[4,estim] = stderror;
    elseif param == 2;  cpsave[1,estim] = theta; cpsave[2,estim] = cmdtest;  cpsave[3,estim] = cdfchic(cmdtest,3); cpsave[4,estim] = stderror;
    elseif param == 3;  ctsave[1,estim] = theta; ctsave[2,estim] = cmdtest;  ctsave[3,estim] = cdfchic(cmdtest,3); ctsave[4,estim] = stderror;
    elseif param ==nz+4;casave[1,estim] = theta; casave[2,estim] = cmdtest;  casave[3,estim] = cdfchic(cmdtest,3); casave[4,estim] = stderror;
    endif;

    for qq(1, nz, 1);
    if param == qq+3;
      cdsave[(qq-1)*4+1,estim] = theta;
      cdsave[(qq-1)*4+2,estim] = cmdtest;
      cdsave[(qq-1)*4+3,estim] = cdfchic(cmdtest,3);
      cdsave[(qq-1)*4+4,estim] = stderror;
    endif;
    endfor;

estim = estim + 1;
endo;

param = param + 1;
endo;


param = 1; @ 1 is beta, 2 through nz+1 are the alphas,
             nz+2 is R2, and the end is the intercept @
do while param <= nz+3;

  if param == 1;
    saveit = bsave[1,.];
    w = invit(moment(ibsave,0)/n);
  elseif param == nz+2;
    saveit = rsave[1,.];
    w = invit(moment(irsave,0)/n);
  elseif param == nz+3;
    saveit = isave[1,.];
    w = invit(moment(iisave,0)/n);
  endif;


    for qq(1, nz, 1);
    if param == qq+1;
      saveit = zsave[(2*qq-1),.];
      w = invit(moment(izsave[((qq-1)*n+1):(qq*n),.],0)/n);
    endif;
    endfor;

    g=ones(nyr,1);
    theta = invpd(g'w*g)*g'w*saveit';

@ Compute t-ratios, etc.                                            @

    f=saveit'-theta;
    stderror = sqrt(invpd(g'w*g)/n);
    cmdtest = n*f'w*f;

    if param == 1;
      bbsave[1,1] = theta;
      bbsave[2,1] = cmdtest;
      bbsave[3,1] = cdfchic(cmdtest,3);
      bbsave[4,1] = stderror;
    elseif param == nz+2;
      brsave[1,1] = theta;
      brsave[2,1] = cmdtest;
      brsave[3,1] = cdfchic(cmdtest,3);
      brsave[4,1] = stderror;
    elseif param == nz+3;
      bisave[1,1] = theta;
      bisave[2,1] = cmdtest;
      bisave[3,1] = cdfchic(cmdtest,3);
      bisave[4,1] = stderror;
    endif;

    for qq(1, nz, 1);
    if param == qq+1;
      bzsave[(qq-1)*4+1,1] = theta;
      bzsave[(qq-1)*4+2,1] = cmdtest;
      bzsave[(qq-1)*4+3,1] = cdfchic(cmdtest,3);
      bzsave[(qq-1)*4+4,1] = stderror;
    endif;
    endfor;

param = param + 1;
endo;



@----------------------------------------------------------------------------------------------
This section computes the sum of the alpha_i's and their standard errors for the CMD estimators.
----------------------------------------------------------------------------------------------@


  idx = seqa(1,nestim,nyr) - 1;

estim = 1;
do while estim <= nestim;

  idx = idx+1;
  saveit = sumc(dsave[seqa(1,2,nz-(ndum*strup)),idx]);
  inflsave = zeros(n,nyr);
  for qq(1, nz-(ndum*strup), 1);
  inflsave = inflsave + idsave[((qq-1)*n+1):(qq*n),idx];
  endfor;
  w = invit(moment(inflsave,0)/n);
  g=ones(nyr,1);
  theta = invpd(g'w*g)*g'w*saveit;

@ Save the intermediate results. @

  disave[1,idx] = saveit';
  disave[2,idx]= (sqrt(diag(moment(inflsave,0)./n^2)))';

@ Compute t-ratios, etc.                                            @

    f=saveit-theta;
    stderror = sqrt(invpd(g'w*g)/n);
    cmdtest = n*f'w*f;

      cdisave[1,estim] = theta;
      cdisave[2,estim] = cmdtest;
      cdisave[3,estim] = cdfchic(cmdtest,3);
      cdisave[4,estim] = stderror;

estim = estim + 1;
endo;


  saveit = sumc(zsave[seqa(1,2,nz-(ndum*strup)),.]);
  inflsave = zeros(n,nyr);
  for qq(1, nz-(ndum*strup), 1);
  inflsave = inflsave + izsave[((qq-1)*n+1):(qq*n),.];
  endfor;
  w = invit(moment(inflsave,0)/n);
  g=ones(nyr,1);
  theta = invpd(g'w*g)*g'w*saveit;

@ Save the intermediate results. @

  zisave[1,.] = saveit';
  zisave[2,.]= (sqrt(diag(moment(inflsave,0)./n^2)))';

@ Compute t-ratios, etc.                                            @

    f=saveit-theta;
    stderror = sqrt(invpd(g'w*g)/n);
    cmdtest = n*f'w*f;

      bzisave[1,1] = theta;
      bzisave[2,1] = cmdtest;
      bzisave[3,1] = cdfchic(cmdtest,3);
      bzisave[4,1] = stderror;

@----------------------------------------------------------------------------------------------
This section computes the sum of the straight dummies
and their standard errors for the CMD estimators.
----------------------------------------------------------------------------------------------@


if strup == 1;

  idx = seqa(1,nestim,nyr) - 1;

estim = 1;
do while estim <= nestim;

  idx = idx+1;
  saveit = asave[1,idx]' + sumc(dsave[seqa(2*(nz-(ndum*strup))+1,2,ndum),idx]);
  inflsave = zeros(n,nyr);
  for qq(nz-(ndum*strup)+1, nz, 1);
  inflsave = inflsave + idsave[((qq-1)*n+1):(qq*n),idx];
  endfor;
  inflsave = inflsave + iasave[1:n,idx];

  w = invit(moment(inflsave,0)/n);
  g=ones(nyr,1);
  theta = invpd(g'w*g)*g'w*saveit;

@ Save the intermediate results. @

  ddsave[1,idx] = saveit';
  ddsave[2,idx]= (sqrt(diag(moment(inflsave,0)./n^2)))';

@ Compute t-ratios, etc.                                            @

    f=saveit-theta;
    stderror = sqrt(invpd(g'w*g)/n);
    cmdtest = n*f'w*f;

      cddsave[1,estim] = theta;
      cddsave[2,estim] = cmdtest;
      cddsave[3,estim] = cdfchic(cmdtest,3);
      cddsave[4,estim] = stderror;

estim = estim + 1;
endo;


  saveit = isave[1,.]' + sumc(zsave[seqa(2*(nz-(ndum*strup))+1,2,ndum),.]);
  inflsave = zeros(n,nyr);
  for qq(nz-(ndum*strup)+1, nz, 1);
  inflsave = inflsave + izsave[((qq-1)*n+1):(qq*n),.];
  endfor;
  inflsave = inflsave + iisave[1:n,.];
  w = invit(moment(inflsave,0)/n);
  g=ones(nyr,1);
  theta = invpd(g'w*g)*g'w*saveit;

@ Save the intermediate results. @

  zdsave[1,.] = saveit';
  zdsave[2,.]= (sqrt(diag(moment(inflsave,0)./n^2)))';

@ Compute t-ratios, etc.                                            @

    f=saveit-theta;
    stderror = sqrt(invpd(g'w*g)/n);
    cmdtest = n*f'w*f;

      bzdsave[1,1] = theta;
      bzdsave[2,1] = cmdtest;
      bzdsave[3,1] = cdfchic(cmdtest,3);
      bzdsave[4,1] = stderror;

endif;

format /rd 8,3;

/*  This business makes a LaTeX table.   */

print;
format /rd 8,3;

print;

for jj(1, nyr, 1);
"& " isave[1,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " asave[1,qq];; endfor; " \\" "\\";
"&(" isave[2,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); ")&(" asave[2,qq];; endfor; ") \\" "\\";
endfor;
"& " bisave[1,1];; for qq(1, nestim, 1); " & " casave[1,qq];; endfor; " \\" "\\";
"&(" bisave[4,1];; for qq(1, nestim, 1); ")&(" casave[4,qq];; endfor; ") \\" "\\";
"&[" bisave[1,1]/bisave[4,1];; for qq(1, nestim, 1); "]&[" casave[1,qq]/casave[4,qq];; endfor; "] \\" "\\";
"& " bisave[2,1];; for qq(1, nestim, 1); " & " casave[2,qq];; endfor; " \\" "\\";
"& " bisave[3,1];; for qq(1, nestim, 1); " & " casave[3,qq];; endfor; " \\" "\\";
"& " meanc(isave[1,.]');; idx=seqa(1,nestim,nyr); for qq(1, nestim, 1); " & " meanc(asave[1,idx+qq-1]');; endfor; " \\" "\\";
print; print;


for jj(1, nyr, 1);
"& " bsave[1,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " csave[1,qq];; endfor; " \\" "\\";
"&(" bsave[2,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); ")&(" csave[2,qq];; endfor; ") \\" "\\";
endfor;
"& " bbsave[1,1];; for qq(1, nestim, 1); " & " ccsave[1,qq];; endfor; " \\" "\\";
"&(" bbsave[4,1];; for qq(1, nestim, 1); ")&(" ccsave[4,qq];; endfor; ") \\" "\\";
"&[" bbsave[1,1]/bbsave[4,1];; for qq(1, nestim, 1); "]&[" ccsave[1,qq]/ccsave[4,qq];; endfor; "] \\" "\\";
"& " bbsave[2,1];; for qq(1, nestim, 1); " & " ccsave[2,qq];; endfor; " \\" "\\";
"& " bbsave[3,1];; for qq(1, nestim, 1); " & " ccsave[3,qq];; endfor; " \\" "\\";
"& " meanc(bsave[1,.]');; idx=seqa(1,nestim,nyr); for qq(1, nestim, 1); " & " meanc(csave[1,idx+qq-1]');; endfor; " \\" "\\";
print; print;

for ee(1, nz, 1);

for jj(1, nyr, 1);
"& " zsave[(ee-1)*2+1,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " dsave[(ee-1)*2+1,qq];; endfor; " \\" "\\";
"&(" zsave[ee*2,jj];;       for qq(((jj-1)*nestim+1), (nestim*jj), 1); ")&(" dsave[ee*2,qq];; endfor; ") \\" "\\";
endfor;
"& " bzsave[(ee-1)*4+1,1];; for qq(1, nestim, 1); " & " cdsave[(ee-1)*4+1,qq];; endfor; " \\" "\\";
"&(" bzsave[(ee-1)*4+4,1];; for qq(1, nestim, 1); ")&(" cdsave[(ee-1)*4+4,qq];; endfor; ") \\" "\\";
"&[" bzsave[(ee-1)*4+1,1]/bzsave[(ee-1)*4+4,1];; for qq(1, nestim, 1); "]&[" cdsave[(ee-1)*4+1,qq]/cdsave[(ee-1)*4+4,qq];; endfor; "] \\" "\\";
"& " bzsave[(ee-1)*4+2,1];; for qq(1, nestim, 1); " & " cdsave[(ee-1)*4+2,qq];; endfor; " \\" "\\";
"& " bzsave[(ee-1)*4+3,1];; for qq(1, nestim, 1); " & " cdsave[(ee-1)*4+3,qq];; endfor; " \\" "\\";
"& " meanc(zsave[(ee-1)*2+1,.]');; idx=seqa(1,nestim,nyr); for qq(1, nestim, 1); " & " meanc(dsave[(ee-1)*2+1,idx+qq-1]');; endfor; " \\" "\\";
print; print;

endfor;                                      print;print;
  "Sum of the interaction terms";print;print;
for jj(1, nyr, 1);
"& " zisave[1,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " disave[1,qq];; endfor; " \\" "\\";
"&(" zisave[2,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); ")&(" disave[2,qq];; endfor; ") \\" "\\";
endfor;
"& " bzisave[1,1];; for qq(1, nestim, 1); " & " cdisave[1,qq];; endfor; " \\" "\\";
"&(" bzisave[4,1];; for qq(1, nestim, 1); ")&(" cdisave[4,qq];; endfor; ") \\" "\\";
"&[" bzisave[1,1]/bzisave[4,1];; for qq(1, nestim, 1); "]&[" cdisave[1,qq]/cdisave[4,qq];; endfor; "] \\" "\\";
"& " bzisave[2,1];; for qq(1, nestim, 1); " & " cdisave[2,qq];; endfor; " \\" "\\";
"& " bzisave[3,1];; for qq(1, nestim, 1); " & " cdisave[3,qq];; endfor; " \\" "\\";
"& " meanc(zisave[1,.]');; idx=seqa(1,nestim,nyr); for qq(1, nestim, 1); " & " meanc(disave[1,idx+qq-1]');; endfor; " \\" "\\";
print; print;

if strup == 1;                                     print;print;
                  "sum of the dummies";print;print;
for jj(1, nyr, 1);
"& " zdsave[1,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " ddsave[1,qq];; endfor; " \\" "\\";
"&(" zdsave[2,jj];; for qq(((jj-1)*nestim+1), (nestim*jj), 1); ")&(" ddsave[2,qq];; endfor; ") \\" "\\";
endfor;
"& " bzdsave[1,1];; for qq(1, nestim, 1); " & " cddsave[1,qq];; endfor; " \\" "\\";
"&(" bzdsave[4,1];; for qq(1, nestim, 1); ")&(" cddsave[4,qq];; endfor; ") \\" "\\";
"&[" bzdsave[1,1]/bzdsave[4,1];; for qq(1, nestim, 1); "]&[" cddsave[1,qq]/cddsave[4,qq];; endfor; "] \\" "\\";
"& " bzdsave[2,1];; for qq(1, nestim, 1); " & " cddsave[2,qq];; endfor; " \\" "\\";
"& " bzdsave[3,1];; for qq(1, nestim, 1); " & " cddsave[3,qq];; endfor; " \\" "\\";
"& " meanc(zdsave[1,.]');; idx=seqa(1,nestim,nyr); for qq(1, nestim, 1); " & " meanc(ddsave[1,idx+qq-1]');; endfor; " \\" "\\";
print; print;

endif;

for jj(1, nyr, 1);
"& " rsave[1,jj];;      for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " psave[1,qq];; endfor; " \\" "\\";
"&(" rsave[2,jj];;")";; for qq(((jj-1)*nestim+1), (nestim*jj), 1); "&(" psave[2,qq];; ")";; endfor; " \\" "\\";
endfor;
"& " brsave[1,1];;   for qq(1, nestim, 1); " & " cpsave[1,qq];; endfor; " \\" "\\";
"&(" brsave[4,1];;")";;  for qq(1, nestim, 1); "&(" cpsave[4,qq];; ")";;  endfor; " \\" "\\";
"&[" brsave[1,1]/brsave[4,1];; for qq(1, nestim, 1); "]&[" cpsave[1,qq]/cpsave[4,qq];; endfor; "] \\" "\\";
"& " brsave[2,1];;  for qq(1, nestim, 1); " & " cpsave[2,qq];; endfor; " \\" "\\";
"& " brsave[3,1];;  for qq(1, nestim, 1); " & " cpsave[3,qq];; endfor; " \\" "\\";
"& " meanc(rsave[1,.]');;   idx=seqa(1,nestim,nyr); for qq(1, nestim, 1); " & " meanc(psave[1,idx+qq-1]');; endfor; " \\" "\\";
print; print;

for jj(1, nyr, 1);
for qq(((jj-1)*nestim+1), (nestim*jj), 1); " & " tsave[1,qq];; endfor; " \\" "\\";
" ";;for qq(((jj-1)*nestim+1), (nestim*jj), 1); "&(" tsave[2,qq];; ")";; endfor; " \\" "\\";
endfor;
for qq(1, nestim, 1); " & " ctsave[1,qq];; endfor; " \\" "\\";
" ";;for qq(1, nestim, 1); "&(" ctsave[4,qq];; ")";;  endfor; " \\" "\\";
" ";;for qq(1, nestim, 1); "&[" ctsave[1,qq]/ctsave[4,qq];; "]";;  endfor; " \\" "\\";
for qq(1, nestim, 1); " & " ctsave[2,qq];; endfor; " \\" "\\";
for qq(1, nestim, 1); " & " ctsave[3,qq];; endfor; " \\" "\\";
idx=seqa(1,nestim,nyr); for qq(1, nestim, 1); " & " meanc(tsave[1,idx+qq-1]');; endfor; " \\" "\\";
print; print;


for jj(1, nyr, 1);
for qq(((jj-1)*nestim+2), (nestim*jj), 1); " & " chisave[1,qq];; endfor; " \\" "\\";
endfor;
print;print;


for jj(1, nyr, 1);
for qq(((jj-1)*nestim+2), (nestim*jj), 1); " & " chisave[2,qq];; endfor; " \\" "\\";
endfor;
print;print;


/*----------------This business makes another LaTeX table.---------------------*/


"& " bisave[1,1];; for qq(1, nestim, 1); " & " casave[1,qq];; endfor; " & " bzdsave[1,1];;  for qq(1, nestim, 1); " & " cddsave[1,qq];; endfor; " \\" "\\";
"&(" bisave[4,1];; for qq(1, nestim, 1); ")&(" casave[4,qq];; endfor; ")&(" bzdsave[4,1];;  for qq(1, nestim, 1); ")&(" cddsave[4,qq];; endfor; ") \\" "\\";
print; print;

"& " bzsave[1,1];; for qq(1, nestim, 1); " & " cdsave[1,qq];; endfor; " & " bzisave[1,1];;  for qq(1, nestim, 1); " & " cdisave[1,qq];; endfor; " \\" "\\";
"&(" bzsave[4,1];; for qq(1, nestim, 1); ")&(" cdsave[4,qq];; endfor; ")&(" bzisave[4,1];;  for qq(1, nestim, 1); ")&(" cdisave[4,qq];; endfor; ") \\" "\\";
print; print;

"& " bbsave[1,1];; for qq(1, nestim, 1); " & " ccsave[1,qq];; endfor;  " \\" "\\";
"&(" bbsave[4,1];; for qq(1, nestim, 1); ")&(" ccsave[4,qq];; endfor;  ") \\" "\\";
print; print;

for qq(1, nestim, 1); " & " cpsave[1,qq];; endfor; for qq(1, nestim, 1); " & " ctsave[1,qq];; endfor; " \\" "\\";
for qq(1, nestim, 1); ")&(" cpsave[4,qq];; endfor; for qq(1, nestim, 1); ")&(" ctsave[4,qq];; endfor; ") \\" "\\";
print; print;

"-------------------------------------------------------------------------------";
"-------------------------------------------------------------------------------";
"-------------------------------------------------------------------------------";
"-------------------------------------------------------------------------------";
"-------------------------------------------------------------------------------";

pick = pick + 1;
endo;

ndum = ndum - 1;
endo;

strup = strup + 1;
endo;


goto eop;

@ ----------------- Subroutines follow --------------- @

@ Subroutine to define the f vector for the optimal weighting matrix. @

proc optw(a,mmm,ooo,wf);
local f, ei;

f = zeros(n,neq);

f[.,1] = varget(mmm[1]) - varget(ooo[1]);
f[.,2] = varget(mmm[2]) - varget(ooo[2]);
f[.,3] = varget(mmm[3]) - varget(ooo[3]);
f[.,4] = varget(mmm[4]) - varget(ooo[4]);
f[.,5] = varget(mmm[5]) - varget(ooo[5]);
if estim > 1;
f[.,6] = varget(mmm[6]) - varget(ooo[6]);
f[.,7] = varget(mmm[7]) - varget(ooo[7]);
f[.,8] = varget(mmm[8]) - varget(ooo[8]);
endif;
if estim > 2;
f[.,9] = varget(mmm[9]) - varget(ooo[9]);
f[.,10] = varget(mmm[10]) - varget(ooo[10]);
f[.,11] = varget(mmm[11]) - varget(ooo[11]);
f[.,12] = varget(mmm[12]) - varget(ooo[12]);
f[.,13] = varget(mmm[13]) - varget(ooo[13]);
f[.,14] = varget(mmm[14]) - varget(ooo[14]);
endif;
if estim > 3;
f[.,15] = varget(mmm[15]) - varget(ooo[15]);
f[.,16] = varget(mmm[16]) - varget(ooo[16]);
f[.,17] = varget(mmm[17]) - varget(ooo[17]);
f[.,18] = varget(mmm[18]) - varget(ooo[18]);
f[.,19] = varget(mmm[19]) - varget(ooo[19]);
f[.,20] = varget(mmm[20]) - varget(ooo[20]);
f[.,21] = varget(mmm[21]) - varget(ooo[21]);
endif;
if estim > 4;
f[.,22] = varget(mmm[22]) - varget(ooo[22]);
f[.,23] = varget(mmm[23]) - varget(ooo[23]);
f[.,24] = varget(mmm[24]) - varget(ooo[24]);
f[.,25] = varget(mmm[25]) - varget(ooo[25]);
f[.,26] = varget(mmm[26]) - varget(ooo[26]);
f[.,27] = varget(mmm[27]) - varget(ooo[27]);
f[.,28] = varget(mmm[28]) - varget(ooo[28]);
f[.,29] = varget(mmm[29]) - varget(ooo[29]);
endif;
@-------------------------------------------------------@
@  This part makes the standard error adjustment.       @
@-------------------------------------------------------@


if wf == 2;

ei = zeros(n,neq);

ei[.,4] = (-2*Eyx_z'inEzz*zy_d' - Ey2_z'inEzz*zx_d')';
ei[.,5] = (-Ex2_z'inEzz*zy_d' - 2*Eyx_z'inEzz*zx_d')';

if estim > 1;
ei[.,6] = (-3*Ey2x_z'inEzz*zy_d' - Ey3_z'inEzz*zx_d')';
ei[.,7] = (-2*Eyx2_z'inEzz*zy_d' - 2*Ey2x_z'inEzz*zx_d')';
ei[.,8] = (-Ex3_z'inEzz*zy_d' - 3*Eyx2_z'inEzz*zx_d')';
endif;

if estim > 2;
ei[.,9]  = (-3*Ex2_z'inEzz*zx_d')';
ei[.,10] = (-3*Ey2_z'inEzz*zy_d')';
ei[.,11] = (-4*Ey3x_z'inEzz*zy_d' - Ey4_z'inEzz*zx_d')';
ei[.,12] = (-3*Ey2x2_z'inEzz*zy_d' - 2*Ey3x_z'inEzz*zx_d')';
ei[.,13] = (-2*Eyx3_z'inEzz*zy_d' - 3*Ey2x2_z'inEzz*zx_d')';
ei[.,14] = (-Ex4_z'inEzz*zy_d' - 4*Eyx3_z'inEzz*zx_d')';
endif;

if estim > 3;
ei[.,15] = (-4*Ey3_z'inEzz*zy_d')';
ei[.,16] = (-4*Ex3_z'inEzz*zx_d')';

ei[.,17] = (-5*Ey4x_z'inEzz*zy_d' -  Ey5_z'inEzz*zx_d')';
ei[.,18] = (-4*Ey3x2_z'inEzz*zy_d' - 2*Ey4x_z'inEzz*zx_d')';
ei[.,19] = (-3*Ey2x3_z'inEzz*zy_d' - 3*Ey3x2_z'inEzz*zx_d')';
ei[.,20] = (-2*Eyx4_z'inEzz*zy_d' - 4*Ey2x3_z'inEzz*zx_d')';
ei[.,21] = (-Ex5_z'inEzz*zy_d' -   5*Eyx4_z'inEzz*zx_d')';
endif;

if estim > 4;
ei[.,22] = (-5*Ey4_z'inEzz*zy_d')';
ei[.,23] = (-5*Ex4_z'inEzz*zx_d')';
ei[.,24] = (-6*Ey5x_z'inEzz*zy_d' -     Ey6_z'inEzz*zx_d')';
ei[.,25] = (-5*Ey4x2_z'inEzz*zy_d' -  2*Ey5x_z'inEzz*zx_d')';
ei[.,26] = (-4*Ey3x3_z'inEzz*zy_d' -  3*Ey4x2_z'inEzz*zx_d')';
ei[.,27] = (-3*Ey2x4_z'inEzz*zy_d' -  4*Ey3x3_z'inEzz*zx_d')';
ei[.,28] = (-2*Eyx5_z'inEzz*zy_d' -   5*Ey2x4_z'inEzz*zx_d')';
ei[.,29] = (-Ex6_z'inEzz*zy_d'   -    6*Eyx5_z'inEzz*zx_d')';
endif;

f = f + ei;
endif;

retp(f);
clear f;
endp;

@ Subroutine to define the f vector. @

proc deff(a,mmm,wf);
local f, ei;

f = zeros(1,neq);

f[.,1] = varget(mmm[1]) - (a[1,1]^2)*a[2,1] - a[3,1];
f[.,2] = varget(mmm[2]) - a[1,1]*a[2,1];
f[.,3] = varget(mmm[3]) - a[2,1] - a[4,1];
f[.,4] = varget(mmm[4]) - (a[1,1]^2)*a[5,1];
f[.,5] = varget(mmm[5]) - a[1,1]*a[5,1];

if estim > 1;

f[.,6] = varget(mmm[6]) - (a[1,1]^3)*a[6,1] - 3*a[1,1]*a[2,1]*a[3,1];
f[.,7] = varget(mmm[7]) - (a[1,1]^2)*( a[6,1] + a[2,1]*a[4,1] )
               - a[3,1]*( a[2,1] + a[4,1] );
f[.,8] = varget(mmm[8]) - a[1,1]*a[6,1] - 3*a[1,1]*a[2,1]*a[4,1];

endif;
if estim > 2;

f[.,9] = varget(mmm[9]) - a[5,1] - a[7,1];

f[.,10] = varget(mmm[10]) - (a[1,1]^3)*a[5,1] - a[8,1];

f[.,11] = varget(mmm[11]) - (a[1,1]^4)*a[9,1] - 6*(a[1,1]^2)*a[5,1]*a[3,1]
               - 4*a[1,1]*a[2,1]*a[8,1];

f[.,12] = varget(mmm[12]) - (a[1,1]^3)*( a[9,1] + a[5,1]*a[4,1] )
                - 3*a[1,1]*a[5,1]*a[3,1] - a[8,1]*( a[2,1] + a[4,1] );

f[.,13] = varget(mmm[13]) - (a[1,1]^2)*(a[9,1] +3*a[5,1]*a[4,1] + a[2,1]*a[7,1])
                - a[3,1]*( a[5,1] + a[7,1] );

f[.,14] = varget(mmm[14]) - a[1,1]*(a[9,1] + 6*a[5,1]*a[4,1] + 4*a[2,1]*a[7,1]);

endif;
if estim > 3;

f[.,15] = varget(mmm[15]) - (a[1,1]^4)*a[6,1] - 6*(a[1,1]^2)*a[2,1]*a[3,1]
                          - a[10,1];

f[.,16] = varget(mmm[16]) - a[6,1] - 6*a[2,1]*a[4,1] - a[11,1];

f[.,17] = varget(mmm[17]) - (a[1,1]^5)*a[12,1] - 10*(a[1,1]^3)*a[6,1]*a[3,1]
                     - 10*(a[1,1]^2)*a[5,1]*a[8,1] - 5*a[1,1]*a[2,1]*a[10,1];

f[.,18] = varget(mmm[18]) - (a[1,1]^4)*( a[12,1] + a[6,1]*a[4,1] )
                       - 6*(a[1,1]^2)*( a[6,1]*a[3,1] + a[2,1]*a[3,1]*a[4,1] )
                       - 4*a[1,1]*a[5,1]*a[8,1] - a[10,1]*( a[2,1] + a[4,1] );

f[.,19] = varget(mmm[19]) - (a[1,1]^3)*(a[12,1]+3*a[6,1]*a[4,1]+ a[5,1]*a[7,1])
                          - a[1,1]*( 3*a[6,1]*a[3,1] + 9*a[2,1]*a[3,1]*a[4,1] )
                          - a[8,1]*( a[5,1] + a[7,1] );

f[.,20] = varget(mmm[20]) - (a[1,1]^2)*( a[12,1] + 6*a[6,1]*a[4,1]
                             + 4*a[5,1]*a[7,1] + a[2,1]*a[11,1] )
                           - a[3,1]*( a[6,1] + 6*a[2,1]*a[4,1] + a[11,1] );

f[.,21] = varget(mmm[21]) - a[1,1]*( a[12,1] + 10*a[6,1]*a[4,1]
                                     + 10*a[5,1]*a[7,1] + 5*a[2,1]*a[11,1] );

endif;
if estim > 4;

f[.,22] = varget(mmm[22]) - (a[1,1]^5)*a[9,1] - 10*(a[1,1]^3)*a[5,1]*a[3,1]
              - 10*(a[1,1]^2)*a[2,1]*a[8,1] - a[13,1];


f[.,23] = varget(mmm[23])
         - a[9,1] - 10*a[5,1]*a[4,1] - 10*a[2,1]*a[7,1] - a[14,1];

f[.,24] = varget(mmm[24]) - (a[1,1]^6)*a[15,1]
               - 15*(a[1,1]^4)*a[9,1]*a[3,1]
               - 20*(a[1,1]^3)*a[6,1]*a[8,1]
               - 15*(a[1,1]^2)*a[5,1]*a[10,1]
               - 6*a[1,1]*a[2,1]*a[13,1];

f[.,25] = varget(mmm[25]) - (a[1,1]^5)*( a[15,1] + a[9,1]*a[4,1] )
                - 10*(a[1,1]^3)*( a[9,1]*a[3,1] + a[5,1]*a[3,1]*a[4,1] )
                - 10*(a[1,1]^2)*( a[6,1]*a[8,1] + a[2,1]*a[8,1]*a[4,1] )
                - 5*a[1,1]*a[5,1]*a[10,1]
                - a[13,1]*( a[2,1] + a[4,1] );


f[.,26] = varget(mmm[26])
- (a[1,1]^4)*( a[15,1] + 3*a[9,1]*a[4,1] + a[6,1]*a[7,1] )
- 6*(a[1,1]^2)*( a[9,1]*a[3,1] + 3*a[5,1]*a[3,1]*a[4,1] + a[2,1]*a[3,1]*a[7,1] )
- 4*a[1,1]*( a[6,1]*a[8,1] + 3*a[2,1]*a[8,1]*a[4,1] )
- a[10,1]*( a[5,1] + a[7,1] );


f[.,27] = varget(mmm[27])
- (a[1,1]^3)*( a[15,1] + 6*a[9,1]*a[4,1] + 4*a[6,1]*a[7,1] + a[5,1]*a[11,1] )
- a[1,1]*( 3*a[9,1]*a[3,1] + 18*a[5,1]*a[3,1]*a[4,1] + 12*a[2,1]*a[3,1]*a[7,1] )
-a[8,1]*( a[6,1] + 6*a[2,1]*a[4,1] + a[11,1] );

f[.,28] = varget(mmm[28])
- (a[1,1]^2)*( a[15,1] + 10*a[9,1]*a[4,1] + 10*a[6,1]*a[7,1] + 5*a[5,1]*a[11,1]
+ a[2,1]*a[14,1] )
- a[3,1]*( a[9,1] + 10*a[5,1]*a[4,1] + 10*a[2,1]*a[7,1] + a[14,1] );

f[.,29] = varget(mmm[29])
- a[1,1]*( a[15,1] + 15*a[9,1]*a[4,1] + 20*a[6,1]*a[7,1] +15*a[5,1]*a[11,1]
+ 6*a[2,1]*a[14,1] );

endif;

if wf == 0; f = f'; endif;

retp(f);
clear f;
endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad1(a);

local dvdc1,dvdc2,dvdc3,dvdc4,dvdc5;
dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];

dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;

dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;

dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;

dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];


retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5;
endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad2(a);

local dvdc1,dvdc2,dvdc3,dvdc4,dvdc5,dvdc6;
dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);
dvdc6  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];
dvdc1[6,1] = -3*(a[1,1]^2)*a[6,1] - 3*a[2,1]*a[3,1];
dvdc1[7,1] = -2*a[1,1]*( a[6,1] + a[2,1]*a[4,1] );
dvdc1[8,1] = -a[6,1] - 3*a[2,1]*a[4,1];

dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;
dvdc2[6,1] = -3*a[1,1]*a[3,1];
dvdc2[7,1] = -(a[1,1]^2)*a[4,1] - a[3,1];
dvdc2[8,1] = -3*a[1,1]*a[4,1];

dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;
dvdc3[6,1] = -3*a[1,1]*a[2,1];
dvdc3[7,1] = -( a[2,1] + a[4,1] );
dvdc3[8,1] = 0;

dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;
dvdc4[6,1] = 0;
dvdc4[7,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc4[8,1] = -3*a[1,1]*a[2,1];

dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];
dvdc5[6,1] = 0;
dvdc5[7,1] = 0;
dvdc5[8,1] = 0;

dvdc6[1,1] = 0;
dvdc6[2,1] = 0;
dvdc6[3,1] = 0;
dvdc6[4,1] = 0;
dvdc6[5,1] = 0;
dvdc6[6,1] = -(a[1,1]^3);
dvdc6[7,1] = -(a[1,1]^2);
dvdc6[8,1] = -a[1,1];


retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5~dvdc6);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5, dvdc6;
endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad3(a);

local dvdc,dvdc1,dvdc2,dvdc3,dvdc4,dvdc5,dvdc6,dvdc7,dvdc8,dvdc9;

dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);
dvdc6  = zeros(neq,1);
dvdc7  = zeros(neq,1);
dvdc8  = zeros(neq,1);
dvdc9  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];
dvdc1[6,1] = -3*(a[1,1]^2)*a[6,1] - 3*a[2,1]*a[3,1];
dvdc1[7,1] = -2*a[1,1]*( a[6,1] + a[2,1]*a[4,1] );
dvdc1[8,1] = -a[6,1] - 3*a[2,1]*a[4,1];
dvdc1[9,1] = 0;
dvdc1[10,1] = -3*(a[1,1]^2)*a[5,1];
dvdc1[11,1] = -4*(a[1,1]^3)*a[9,1] - 12*a[1,1]*a[5,1]*a[3,1]
              - 4*a[2,1]*a[8,1];
dvdc1[12,1] = -3*(a[1,1]^2)*( a[9,1] + a[5,1]*a[4,1] )
              - 3*a[5,1]*a[3,1];
dvdc1[13,1] = -2*a[1,1]*( a[9,1] + 3*a[5,1]*a[4,1] + a[2,1]*a[7,1] );
dvdc1[14,1] = -( a[9,1] + 6*a[5,1]*a[4,1] + 4*a[2,1]*a[7,1] );

dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;
dvdc2[6,1] = -3*a[1,1]*a[3,1];
dvdc2[7,1] = -(a[1,1]^2)*a[4,1] - a[3,1];
dvdc2[8,1] = -3*a[1,1]*a[4,1];
dvdc2[9,1] = 0;
dvdc2[10,1] = 0;
dvdc2[11,1] = -4*a[1,1]*a[8,1];
dvdc2[12,1] = -a[8,1];
dvdc2[13,1] = -(a[1,1]^2)*a[7,1];
dvdc2[14,1] = -4*a[1,1]*a[7,1];

dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;
dvdc3[6,1] = -3*a[1,1]*a[2,1];
dvdc3[7,1] = -( a[2,1] + a[4,1] );
dvdc3[8,1] = 0;
dvdc3[9,1] = 0;
dvdc3[10,1] = 0;
dvdc3[11,1] = -6*(a[1,1]^2)*a[5,1];
dvdc3[12,1] = -3*a[1,1]*a[5,1];
dvdc3[13,1] = -( a[5,1] + a[7,1] );
dvdc3[14,1] = 0;

dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;
dvdc4[6,1] = 0;
dvdc4[7,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc4[8,1] = -3*a[1,1]*a[2,1];
dvdc4[9,1] = 0;
dvdc4[10,1] = 0;
dvdc4[11,1] = 0;
dvdc4[12,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc4[13,1] = -3*(a[1,1]^2)*a[5,1];
dvdc4[14,1] = -6*a[1,1]*a[5,1];

dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];
dvdc5[6,1] = 0;
dvdc5[7,1] = 0;
dvdc5[8,1] = 0;
dvdc5[9,1] = -1;
dvdc5[10,1] = -(a[1,1]^3);
dvdc5[11,1] = -6*(a[1,1]^2)*a[3,1];
dvdc5[12,1] = -(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc5[13,1] = -3*(a[1,1]^2)*a[4,1] -a[3,1];
dvdc5[14,1] = -6*a[1,1]*a[4,1];

dvdc6[1,1] = 0;
dvdc6[2,1] = 0;
dvdc6[3,1] = 0;
dvdc6[4,1] = 0;
dvdc6[5,1] = 0;
dvdc6[6,1] = -(a[1,1]^3);
dvdc6[7,1] = -(a[1,1]^2);
dvdc6[8,1] = -a[1,1];
dvdc6[9,1] = 0;
dvdc6[10,1] = 0;
dvdc6[11,1] = 0;
dvdc6[12,1] = 0;
dvdc6[13,1] = 0;
dvdc6[14,1] = 0;

dvdc7[1,1] = 0;
dvdc7[2,1] = 0;
dvdc7[3,1] = 0;
dvdc7[4,1] = 0;
dvdc7[5,1] = 0;
dvdc7[6,1] = 0;
dvdc7[7,1] = 0;
dvdc7[8,1] = 0;
dvdc7[9,1] = -1;
dvdc7[10,1] = 0;
dvdc7[11,1] = 0;
dvdc7[12,1] = 0;
dvdc7[13,1] = -(a[1,1]^2)*a[2,1] -a[3,1];
dvdc7[14,1] = -4*a[1,1]*a[2,1];

dvdc8[1,1] = 0;
dvdc8[2,1] = 0;
dvdc8[3,1] = 0;
dvdc8[4,1] = 0;
dvdc8[5,1] = 0;
dvdc8[6,1] = 0;
dvdc8[7,1] = 0;
dvdc8[8,1] = 0;
dvdc8[9,1] = 0;
dvdc8[10,1] = -1;
dvdc8[11,1] = -4*a[1,1]*a[2,1];
dvdc8[12,1] = -( a[2,1] + a[4,1] );
dvdc8[13,1] = 0;
dvdc8[14,1] = 0;

dvdc9[1,1] = 0;
dvdc9[2,1] = 0;
dvdc9[3,1] = 0;
dvdc9[4,1] = 0;
dvdc9[5,1] = 0;
dvdc9[6,1] = 0;
dvdc9[7,1] = 0;
dvdc9[8,1] = 0;
dvdc9[9,1] = 0;
dvdc9[10,1] = 0;
dvdc9[11,1] = -(a[1,1]^4);
dvdc9[12,1] = -(a[1,1]^3);
dvdc9[13,1] = -(a[1,1]^2);
dvdc9[14,1] = -a[1,1];


retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5~dvdc6~dvdc7~dvdc8~dvdc9);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5, dvdc6, dvdc7, dvdc8, dvdc9;

endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad4(a);

local dvdc,dvdc1,dvdc2,dvdc3,dvdc4,dvdc5,dvdc6,dvdc7,dvdc8,dvdc9,
dvdc10,dvdc11,dvdc12;

dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);
dvdc6  = zeros(neq,1);
dvdc7  = zeros(neq,1);
dvdc8  = zeros(neq,1);
dvdc9  = zeros(neq,1);
dvdc10  = zeros(neq,1);
dvdc11  = zeros(neq,1);
dvdc12  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];
dvdc1[6,1] = -3*(a[1,1]^2)*a[6,1] - 3*a[2,1]*a[3,1];
dvdc1[7,1] = -2*a[1,1]*( a[6,1] + a[2,1]*a[4,1] );
dvdc1[8,1] = -a[6,1] - 3*a[2,1]*a[4,1];
dvdc1[9,1] = 0;
dvdc1[10,1] = -3*(a[1,1]^2)*a[5,1];
dvdc1[11,1] = -4*(a[1,1]^3)*a[9,1] - 12*a[1,1]*a[5,1]*a[3,1]
              - 4*a[2,1]*a[8,1];
dvdc1[12,1] = -3*(a[1,1]^2)*( a[9,1] + a[5,1]*a[4,1] )
              - 3*a[5,1]*a[3,1];
dvdc1[13,1] = -2*a[1,1]*( a[9,1] + 3*a[5,1]*a[4,1] + a[2,1]*a[7,1] );
dvdc1[14,1] = -( a[9,1] + 6*a[5,1]*a[4,1] + 4*a[2,1]*a[7,1] );

dvdc1[15,1] = -4*(a[1,1]^3)*a[6,1] -12*a[1,1]*a[2,1]*a[3,1];
dvdc1[16,1] = 0;
dvdc1[17,1] = -5*(a[1,1]^4)*a[12,1] - 30*(a[1,1]^2)*a[6,1]*a[3,1]
              -20*a[1,1]*a[5,1]*a[8,1] - 5*a[2,1]*a[10,1];
dvdc1[18,1] = -4*(a[1,1]^3)*( a[12,1] + a[6,1]*a[4,1] )
              -12*a[1,1]*( a[6,1]*a[3,1] + a[2,1]*a[3,1]*a[4,1] )
              -4*a[5,1]*a[8,1];
dvdc1[19,1] = -3*(a[1,1]^2)*( a[12,1] + 3*a[6,1]*a[4,1] + a[5,1]*a[7,1] )
              -( 3*a[6,1]*a[3,1] + 9*a[2,1]*a[3,1]*a[4,1] );
dvdc1[20,1] = -2*a[1,1]*( a[12,1] + 6*a[6,1]*a[4,1] + 4*a[5,1]*a[7,1]
                         + a[2,1]*a[11,1] );
dvdc1[21,1] = - ( a[12,1] + 10*a[6,1]*a[4,1] + 10*a[5,1]*a[7,1]
                  + 5*a[2,1]*a[11,1] );


dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;
dvdc2[6,1] = -3*a[1,1]*a[3,1];
dvdc2[7,1] = -(a[1,1]^2)*a[4,1] - a[3,1];
dvdc2[8,1] = -3*a[1,1]*a[4,1];
dvdc2[9,1] = 0;
dvdc2[10,1] = 0;
dvdc2[11,1] = -4*a[1,1]*a[8,1];
dvdc2[12,1] = -a[8,1];
dvdc2[13,1] = -(a[1,1]^2)*a[7,1];
dvdc2[14,1] = -4*a[1,1]*a[7,1];

dvdc2[15,1] = -6*(a[1,1]^2)*a[3,1];
dvdc2[16,1] = -6*a[4,1];
dvdc2[17,1] = -5*a[1,1]*a[10,1];
dvdc2[18,1] = -6*(a[1,1]^2)*a[3,1]*a[4,1] - a[10,1];
dvdc2[19,1] = -9*a[1,1]*a[3,1]*a[4,1];
dvdc2[20,1] = -(a[1,1]^2)*a[11,1] - 6*a[3,1]*a[4,1];
dvdc2[21,1] = -5*a[1,1]*a[11,1];



dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;
dvdc3[6,1] = -3*a[1,1]*a[2,1];
dvdc3[7,1] = -( a[2,1] + a[4,1] );
dvdc3[8,1] = 0;
dvdc3[9,1] = 0;
dvdc3[10,1] = 0;
dvdc3[11,1] = -6*(a[1,1]^2)*a[5,1];
dvdc3[12,1] = -3*a[1,1]*a[5,1];
dvdc3[13,1] = -( a[5,1] + a[7,1] );
dvdc3[14,1] = 0;

dvdc3[15,1] = -6*(a[1,1]^2)*a[2,1];
dvdc3[16,1] = 0;
dvdc3[17,1] = -10*(a[1,1]^3)*a[6,1];
dvdc3[18,1] = -6*(a[1,1]^2)*( a[6,1] + a[2,1] );
dvdc3[19,1] = -a[1,1]*( 3*a[6,1] + 9*a[2,1]*a[4,1] );
dvdc3[20,1] = -( a[6,1] + 6*a[2,1]*a[4,1] + a[11,1] );
dvdc3[21,1] = 0;



dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;
dvdc4[6,1] = 0;
dvdc4[7,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc4[8,1] = -3*a[1,1]*a[2,1];
dvdc4[9,1] = 0;
dvdc4[10,1] = 0;
dvdc4[11,1] = 0;
dvdc4[12,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc4[13,1] = -3*(a[1,1]^2)*a[5,1];
dvdc4[14,1] = -6*a[1,1]*a[5,1];

dvdc4[15,1] = 0;
dvdc4[16,1] = -6*a[2,1];
dvdc4[17,1] = 0;
dvdc4[18,1] = -(a[1,1]^4)*a[6,1] - 6*(a[1,1]^2)*a[2,1]*a[3,1] - a[10,1];
dvdc4[19,1] = -3*(a[1,1]^3)*a[6,1] - 9*a[1,1]*a[2,1]*a[3,1];
dvdc4[20,1] = -6*(a[1,1]^2)*a[6,1] - 6*a[2,1]*a[3,1];
dvdc4[21,1] = -10*a[1,1]*a[6,1];



dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];
dvdc5[6,1] = 0;
dvdc5[7,1] = 0;
dvdc5[8,1] = 0;
dvdc5[9,1] = -1;
dvdc5[10,1] = -(a[1,1]^3);
dvdc5[11,1] = -6*(a[1,1]^2)*a[3,1];
dvdc5[12,1] = -(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc5[13,1] = -3*(a[1,1]^2)*a[4,1] -a[3,1];
dvdc5[14,1] = -6*a[1,1]*a[4,1];

dvdc5[15,1] = 0;
dvdc5[16,1] = 0;
dvdc5[17,1] = -10*(a[1,1]^2)*a[8,1];
dvdc5[18,1] = -4*a[1,1]*a[8,1];
dvdc5[19,1] = -(a[1,1]^3)*a[7,1] - a[8,1];
dvdc5[20,1] = -4*(a[1,1]^2)*a[7,1];
dvdc5[21,1] = -10*a[1,1]*a[7,1];



dvdc6[1,1] = 0;
dvdc6[2,1] = 0;
dvdc6[3,1] = 0;
dvdc6[4,1] = 0;
dvdc6[5,1] = 0;
dvdc6[6,1] = -(a[1,1]^3);
dvdc6[7,1] = -(a[1,1]^2);
dvdc6[8,1] = -a[1,1];
dvdc6[9,1] = 0;
dvdc6[10,1] = 0;
dvdc6[11,1] = 0;
dvdc6[12,1] = 0;
dvdc6[13,1] = 0;
dvdc6[14,1] = 0;

dvdc6[15,1] = -(a[1,1]^4);
dvdc6[16,1] = -1;
dvdc6[17,1] = -10*(a[1,1]^3)*a[3,1];
dvdc6[18,1] = -(a[1,1]^4)*a[4,1] - 6*(a[1,1]^2)*a[3,1];
dvdc6[19,1] = -3*(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc6[20,1] = -6*(a[1,1]^2)*a[4,1] - a[3,1];
dvdc6[21,1] = -10*a[1,1]*a[4,1];



dvdc7[1,1] = 0;
dvdc7[2,1] = 0;
dvdc7[3,1] = 0;
dvdc7[4,1] = 0;
dvdc7[5,1] = 0;
dvdc7[6,1] = 0;
dvdc7[7,1] = 0;
dvdc7[8,1] = 0;
dvdc7[9,1] = -1;
dvdc7[10,1] = 0;
dvdc7[11,1] = 0;
dvdc7[12,1] = 0;
dvdc7[13,1] = -(a[1,1]^2)*a[2,1] -a[3,1];
dvdc7[14,1] = -4*a[1,1]*a[2,1];

dvdc7[15,1] = 0;
dvdc7[16,1] = 0;
dvdc7[17,1] = 0;
dvdc7[18,1] = 0;
dvdc7[19,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc7[20,1] = -4*(a[1,1]^2)*a[5,1];
dvdc7[21,1] = -10*a[1,1]*a[5,1];



dvdc8[1,1] = 0;
dvdc8[2,1] = 0;
dvdc8[3,1] = 0;
dvdc8[4,1] = 0;
dvdc8[5,1] = 0;
dvdc8[6,1] = 0;
dvdc8[7,1] = 0;
dvdc8[8,1] = 0;
dvdc8[9,1] = 0;
dvdc8[10,1] = -1;
dvdc8[11,1] = -4*a[1,1]*a[2,1];
dvdc8[12,1] = -( a[2,1] + a[4,1] );
dvdc8[13,1] = 0;
dvdc8[14,1] = 0;

dvdc8[15,1] = 0;
dvdc8[16,1] = 0;
dvdc8[17,1] = -10*(a[1,1]^2)*a[5,1];
dvdc8[18,1] = -4*a[1,1]*a[5,1];
dvdc8[19,1] = -( a[5,1] + a[7,1] );
dvdc8[20,1] = 0;
dvdc8[21,1] = 0;



dvdc9[1,1] = 0;
dvdc9[2,1] = 0;
dvdc9[3,1] = 0;
dvdc9[4,1] = 0;
dvdc9[5,1] = 0;
dvdc9[6,1] = 0;
dvdc9[7,1] = 0;
dvdc9[8,1] = 0;
dvdc9[9,1] = 0;
dvdc9[10,1] = 0;
dvdc9[11,1] = -(a[1,1]^4);
dvdc9[12,1] = -(a[1,1]^3);
dvdc9[13,1] = -(a[1,1]^2);
dvdc9[14,1] = -a[1,1];

dvdc9[15,1] = 0;
dvdc9[16,1] = 0;
dvdc9[17,1] = 0;
dvdc9[18,1] = 0;
dvdc9[19,1] = 0;
dvdc9[20,1] = 0;
dvdc9[21,1] = 0;


dvdc10[1,1] = 0;
dvdc10[2,1] = 0;
dvdc10[3,1] = 0;
dvdc10[4,1] = 0;
dvdc10[5,1] = 0;
dvdc10[6,1] = 0;
dvdc10[7,1] = 0;
dvdc10[8,1] = 0;
dvdc10[9,1] = 0;
dvdc10[10,1] = 0;
dvdc10[11,1] = 0;
dvdc10[12,1] = 0;
dvdc10[13,1] = 0;
dvdc10[14,1] = 0;

dvdc10[15,1] = -1;
dvdc10[16,1] = 0;
dvdc10[17,1] = -5*a[1,1]*a[2,1];
dvdc10[18,1] = -( a[2,1] + a[4,1] );
dvdc10[19,1] = 0;
dvdc10[20,1] = 0;
dvdc10[21,1] = 0;


dvdc11[1,1] = 0;
dvdc11[2,1] = 0;
dvdc11[3,1] = 0;
dvdc11[4,1] = 0;
dvdc11[5,1] = 0;
dvdc11[6,1] = 0;
dvdc11[7,1] = 0;
dvdc11[8,1] = 0;
dvdc11[9,1] = 0;
dvdc11[10,1] = 0;
dvdc11[11,1] = 0;
dvdc11[12,1] = 0;
dvdc11[13,1] = 0;
dvdc11[14,1] = 0;

dvdc11[15,1] = 0;
dvdc11[16,1] = -1;
dvdc11[17,1] = 0;
dvdc11[18,1] = 0;
dvdc11[19,1] = 0;
dvdc11[20,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc11[21,1] = -5*a[1,1]*a[2,1];


dvdc12[1,1] = 0;
dvdc12[2,1] = 0;
dvdc12[3,1] = 0;
dvdc12[4,1] = 0;
dvdc12[5,1] = 0;
dvdc12[6,1] = 0;
dvdc12[7,1] = 0;
dvdc12[8,1] = 0;
dvdc12[9,1] = 0;
dvdc12[10,1] = 0;
dvdc12[11,1] = 0;
dvdc12[12,1] = 0;
dvdc12[13,1] = 0;
dvdc12[14,1] = 0;

dvdc12[15,1] = 0;
dvdc12[16,1] = 0;
dvdc12[17,1] = -(a[1,1]^5);
dvdc12[18,1] = -(a[1,1]^4);
dvdc12[19,1] = -(a[1,1]^3);
dvdc12[20,1] = -(a[1,1]^2);
dvdc12[21,1] = -a[1,1];



retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5~dvdc6~dvdc7~dvdc8~dvdc9
          ~dvdc10~dvdc11~dvdc12);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5, dvdc6, dvdc7, dvdc8, dvdc9,
      dvdc10, dvdc11, dvdc12;

endp;

@  ------------------------------------------------------- @
@ Subroutine to compute gradient matrix. @

proc grad5(a);

local dvdc,dvdc1,dvdc2,dvdc3,dvdc4,dvdc5,dvdc6,dvdc7,dvdc8,dvdc9,
dvdc10,dvdc11,dvdc12,dvdc13,dvdc14,dvdc15;

dvdc1  = zeros(neq,1);
dvdc2  = zeros(neq,1);
dvdc3  = zeros(neq,1);
dvdc4  = zeros(neq,1);
dvdc5  = zeros(neq,1);
dvdc6  = zeros(neq,1);
dvdc7  = zeros(neq,1);
dvdc8  = zeros(neq,1);
dvdc9  = zeros(neq,1);
dvdc10  = zeros(neq,1);
dvdc11  = zeros(neq,1);
dvdc12  = zeros(neq,1);
dvdc13  = zeros(neq,1);
dvdc14  = zeros(neq,1);
dvdc15  = zeros(neq,1);

dvdc1[1,1] = -2*a[1,1]*a[2,1];
dvdc1[2,1] = -a[2,1];
dvdc1[3,1] = 0;
dvdc1[4,1] = -2*a[1,1]*a[5,1];
dvdc1[5,1] = -a[5,1];
dvdc1[6,1] = -3*(a[1,1]^2)*a[6,1] - 3*a[2,1]*a[3,1];
dvdc1[7,1] = -2*a[1,1]*( a[6,1] + a[2,1]*a[4,1] );
dvdc1[8,1] = -a[6,1] - 3*a[2,1]*a[4,1];
dvdc1[9,1] = 0;
dvdc1[10,1] = -3*(a[1,1]^2)*a[5,1];
dvdc1[11,1] = -4*(a[1,1]^3)*a[9,1] - 12*a[1,1]*a[5,1]*a[3,1]
              - 4*a[2,1]*a[8,1];
dvdc1[12,1] = -3*(a[1,1]^2)*( a[9,1] + a[5,1]*a[4,1] )
              - 3*a[5,1]*a[3,1];
dvdc1[13,1] = -2*a[1,1]*( a[9,1] + 3*a[5,1]*a[4,1] + a[2,1]*a[7,1] );
dvdc1[14,1] = -( a[9,1] + 6*a[5,1]*a[4,1] + 4*a[2,1]*a[7,1] );
dvdc1[15,1] = -4*(a[1,1]^3)*a[6,1] -12*a[1,1]*a[2,1]*a[3,1];
dvdc1[16,1] = 0;
dvdc1[17,1] = -5*(a[1,1]^4)*a[12,1] - 30*(a[1,1]^2)*a[6,1]*a[3,1]
              -20*a[1,1]*a[5,1]*a[8,1] - 5*a[2,1]*a[10,1];
dvdc1[18,1] = -4*(a[1,1]^3)*( a[12,1] + a[6,1]*a[4,1] )
              -12*a[1,1]*( a[6,1]*a[3,1] + a[2,1]*a[3,1]*a[4,1] )
              -4*a[5,1]*a[8,1];
dvdc1[19,1] = -3*(a[1,1]^2)*( a[12,1] + 3*a[6,1]*a[4,1] + a[5,1]*a[7,1] )
              -( 3*a[6,1]*a[3,1] + 9*a[2,1]*a[3,1]*a[4,1] );
dvdc1[20,1] = -2*a[1,1]*( a[12,1] + 6*a[6,1]*a[4,1] + 4*a[5,1]*a[7,1]
                         + a[2,1]*a[11,1] );
dvdc1[21,1] = - ( a[12,1] + 10*a[6,1]*a[4,1] + 10*a[5,1]*a[7,1]
                  + 5*a[2,1]*a[11,1] );

dvdc1[22,1] = -5*(a[1,1]^4)*a[9,1] - 30*(a[1,1]^2)*a[5,1]*a[3,1]
              - 20*a[1,1]*a[2,1]*a[8,1];
dvdc1[23,1] = 0;
dvdc1[24,1] = -6*(a[1,1]^5)*a[15,1] - 60*(a[1,1]^3)*a[9,1]*a[3,1]
              -60*(a[1,1]^2)*a[6,1]*a[8,1] - 30*a[1,1]*a[5,1]*a[10,1]
              -6*a[2,1]*a[13,1];
dvdc1[25,1] = -5*(a[1,1]^4)*(a[15,1] + a[9,1]*a[4,1])
              -30*(a[1,1]^2)*(a[9,1]*a[3,1] + a[5,1]*a[3,1]*a[4,1])
              -20*a[1,1]*(a[6,1]*a[8,1] + a[2,1]*a[8,1]*a[4,1])
              -5*a[5,1]*a[10,1];

dvdc1[26,1] = -4*(a[1,1]^3)*(a[15,1] + 3*a[9,1]*a[4,1] + a[6,1]*a[7,1])
-12*a[1,1]*(a[9,1]*a[3,1] + 3*a[5,1]*a[3,1]*a[4,1] + a[2,1]*a[3,1]*a[7,1])
-4*(a[6,1]*a[8,1] + 3*a[2,1]*a[8,1]*a[4,1]);

dvdc1[27,1] =
-3*(a[1,1]^2)*(a[15,1] + 6*a[9,1]*a[4,1] + 4*a[6,1]*a[7,1] + a[5,1]*a[11,1])
-(3*a[9,1]*a[3,1] + 18*a[5,1]*a[3,1]*a[4,1] + 12*a[2,1]*a[3,1]*a[7,1]);

dvdc1[28,1] =
-2*a[1,1]*(a[15,1] + 10*a[9,1]*a[4,1] + 10*a[6,1]*a[7,1] + 5*a[5,1]*a[11,1]
+ a[2,1]*a[14,1]);

dvdc1[29,1] = -(a[15,1] + 15*a[9,1]*a[4,1] + 20*a[6,1]*a[7,1]
              + 15*a[5,1]*a[11,1] + 6*a[2,1]*a[14,1]);



dvdc2[1,1] = -(a[1,1]^2);
dvdc2[2,1] = -a[1,1];
dvdc2[3,1] = -1;
dvdc2[4,1] = 0;
dvdc2[5,1] = 0;
dvdc2[6,1] = -3*a[1,1]*a[3,1];
dvdc2[7,1] = -(a[1,1]^2)*a[4,1] - a[3,1];
dvdc2[8,1] = -3*a[1,1]*a[4,1];
dvdc2[9,1] = 0;
dvdc2[10,1] = 0;
dvdc2[11,1] = -4*a[1,1]*a[8,1];
dvdc2[12,1] = -a[8,1];
dvdc2[13,1] = -(a[1,1]^2)*a[7,1];
dvdc2[14,1] = -4*a[1,1]*a[7,1];
dvdc2[15,1] = -6*(a[1,1]^2)*a[3,1];
dvdc2[16,1] = -6*a[4,1];
dvdc2[17,1] = -5*a[1,1]*a[10,1];
dvdc2[18,1] = -6*(a[1,1]^2)*a[3,1]*a[4,1] - a[10,1];
dvdc2[19,1] = -9*a[1,1]*a[3,1]*a[4,1];
dvdc2[20,1] = -(a[1,1]^2)*a[11,1] - 6*a[3,1]*a[4,1];
dvdc2[21,1] = -5*a[1,1]*a[11,1];

dvdc2[22,1] = -10*(a[1,1]^2)*a[8,1];
dvdc2[23,1] = -10*a[7,1];
dvdc2[24,1] = -6*a[1,1]*a[13,1];
dvdc2[25,1] = -10*(a[1,1]^2)*a[8,1]*a[4,1] - a[13,1];
dvdc2[26,1] = -6*(a[1,1]^2)*a[3,1]*a[7,1] - 12*a[1,1]*a[8,1]*a[4,1];
dvdc2[27,1] = -12*a[1,1]*a[3,1]*a[7,1] - 6*a[8,1]*a[4,1];
dvdc2[28,1] = -(a[1,1]^2)*a[14,1] - 10*a[3,1]*a[7,1];
dvdc2[29,1] = -6*a[1,1]*a[14,1];



dvdc3[1,1] = -1;
dvdc3[2,1] = 0;
dvdc3[3,1] = 0;
dvdc3[4,1] = 0;
dvdc3[5,1] = 0;
dvdc3[6,1] = -3*a[1,1]*a[2,1];
dvdc3[7,1] = -( a[2,1] + a[4,1] );
dvdc3[8,1] = 0;
dvdc3[9,1] = 0;
dvdc3[10,1] = 0;
dvdc3[11,1] = -6*(a[1,1]^2)*a[5,1];
dvdc3[12,1] = -3*a[1,1]*a[5,1];
dvdc3[13,1] = -( a[5,1] + a[7,1] );
dvdc3[14,1] = 0;
dvdc3[15,1] = -6*(a[1,1]^2)*a[2,1];
dvdc3[16,1] = 0;
dvdc3[17,1] = -10*(a[1,1]^3)*a[6,1];
dvdc3[18,1] = -6*(a[1,1]^2)*( a[6,1] + a[2,1] );
dvdc3[19,1] = -a[1,1]*( 3*a[6,1] + 9*a[2,1]*a[4,1] );
dvdc3[20,1] = -( a[6,1] + 6*a[2,1]*a[4,1] + a[11,1] );
dvdc3[21,1] = 0;

dvdc3[22,1] = -10*(a[1,1]^3)*a[5,1];
dvdc3[23,1] = 0;
dvdc3[24,1] = -15*(a[1,1]^4)*a[9,1];
dvdc3[25,1] = -10*(a[1,1]^3)*(a[9,1] + a[5,1]*a[4,1]);
dvdc3[26,1] = -6*(a[1,1]^2)*(a[9,1] + 3*a[5,1]*a[4,1] + a[2,1]*a[7,1]);
dvdc3[27,1] = -a[1,1]*(3*a[9,1] + 18*a[5,1]*a[4,1] + 12*a[2,1]*a[7,1]);
dvdc3[28,1] = -(a[9,1] + 10*a[5,1]*a[4,1] + 10*a[2,1]*a[7,1] + a[14,1]);
dvdc3[29,1] = 0;



dvdc4[1,1] = 0;
dvdc4[2,1] = 0;
dvdc4[3,1] = -1;
dvdc4[4,1] = 0;
dvdc4[5,1] = 0;
dvdc4[6,1] = 0;
dvdc4[7,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc4[8,1] = -3*a[1,1]*a[2,1];
dvdc4[9,1] = 0;
dvdc4[10,1] = 0;
dvdc4[11,1] = 0;
dvdc4[12,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc4[13,1] = -3*(a[1,1]^2)*a[5,1];
dvdc4[14,1] = -6*a[1,1]*a[5,1];
dvdc4[15,1] = 0;
dvdc4[16,1] = -6*a[2,1];
dvdc4[17,1] = 0;
dvdc4[18,1] = -(a[1,1]^4)*a[6,1] - 6*(a[1,1]^2)*a[2,1]*a[3,1] - a[10,1];
dvdc4[19,1] = -3*(a[1,1]^3)*a[6,1] - 9*a[1,1]*a[2,1]*a[3,1];
dvdc4[20,1] = -6*(a[1,1]^2)*a[6,1] - 6*a[2,1]*a[3,1];
dvdc4[21,1] = -10*a[1,1]*a[6,1];

dvdc4[22,1] = 0;
dvdc4[23,1] = -10*a[5,1];
dvdc4[24,1] = 0;
dvdc4[25,1] = -(a[1,1]^5)*a[9,1] - 10*(a[1,1]^3)*a[5,1]*a[3,1]
              -10*(a[1,1]^2)*a[2,1]*a[8,1] - a[13,1];
dvdc4[26,1] = -3*(a[1,1]^4)*a[9,1] - 18*(a[1,1]^2)*a[5,1]*a[3,1]
              -12*a[1,1]*a[2,1]*a[8,1];
dvdc4[27,1] = -6*(a[1,1]^3)*a[9,1] - 18*a[1,1]*a[5,1]*a[3,1]
              -6*a[8,1]*a[2,1];
dvdc4[28,1] = -10*(a[1,1]^2)*a[9,1] - 10*a[3,1]*a[5,1];
dvdc4[29,1] = -15*a[1,1]*a[9,1];



dvdc5[1,1] = 0;
dvdc5[2,1] = 0;
dvdc5[3,1] = 0;
dvdc5[4,1] = -(a[1,1]^2);
dvdc5[5,1] = -a[1,1];
dvdc5[6,1] = 0;
dvdc5[7,1] = 0;
dvdc5[8,1] = 0;
dvdc5[9,1] = -1;
dvdc5[10,1] = -(a[1,1]^3);
dvdc5[11,1] = -6*(a[1,1]^2)*a[3,1];
dvdc5[12,1] = -(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc5[13,1] = -3*(a[1,1]^2)*a[4,1] -a[3,1];
dvdc5[14,1] = -6*a[1,1]*a[4,1];
dvdc5[15,1] = 0;
dvdc5[16,1] = 0;
dvdc5[17,1] = -10*(a[1,1]^2)*a[8,1];
dvdc5[18,1] = -4*a[1,1]*a[8,1];
dvdc5[19,1] = -(a[1,1]^3)*a[7,1] - a[8,1];
dvdc5[20,1] = -4*(a[1,1]^2)*a[7,1];
dvdc5[21,1] = -10*a[1,1]*a[7,1];

dvdc5[22,1] = -10*(a[1,1]^3)*a[3,1];
dvdc5[23,1] = -10*a[4,1];
dvdc5[24,1] = -15*(a[1,1]^2)*a[10,1];
dvdc5[25,1] = -10*(a[1,1]^3)*a[3,1]*a[4,1] - 5*a[1,1]*a[10,1];
dvdc5[26,1] = -18*(a[1,1]^2)*a[3,1]*a[4,1] - a[10,1];
dvdc5[27,1] = -(a[1,1]^3)*a[11,1] - 18*a[1,1]*a[3,1]*a[4,1];
dvdc5[28,1] = -5*(a[1,1]^2)*a[11,1] - 10*a[3,1]*a[4,1];
dvdc5[29,1] = -15*a[1,1]*a[11,1];



dvdc6[1,1] = 0;
dvdc6[2,1] = 0;
dvdc6[3,1] = 0;
dvdc6[4,1] = 0;
dvdc6[5,1] = 0;
dvdc6[6,1] = -(a[1,1]^3);
dvdc6[7,1] = -(a[1,1]^2);
dvdc6[8,1] = -a[1,1];
dvdc6[9,1] = 0;
dvdc6[10,1] = 0;
dvdc6[11,1] = 0;
dvdc6[12,1] = 0;
dvdc6[13,1] = 0;
dvdc6[14,1] = 0;
dvdc6[15,1] = -(a[1,1]^4);
dvdc6[16,1] = -1;
dvdc6[17,1] = -10*(a[1,1]^3)*a[3,1];
dvdc6[18,1] = -(a[1,1]^4)*a[4,1] - 6*(a[1,1]^2)*a[3,1];
dvdc6[19,1] = -3*(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc6[20,1] = -6*(a[1,1]^2)*a[4,1] - a[3,1];
dvdc6[21,1] = -10*a[1,1]*a[4,1];

dvdc6[22,1] = 0;
dvdc6[23,1] = 0;
dvdc6[24,1] = -20*(a[1,1]^3)*a[8,1];
dvdc6[25,1] = -10*(a[1,1]^2)*a[8,1];
dvdc6[26,1] = -(a[1,1]^4)*a[7,1] - 4*a[1,1]*a[8,1];
dvdc6[27,1] = -4*(a[1,1]^3)*a[7,1] - a[8,1];
dvdc6[28,1] = -10*(a[1,1]^2)*a[7,1];
dvdc6[29,1] = -20*a[1,1]*a[7,1];



dvdc7[1,1] = 0;
dvdc7[2,1] = 0;
dvdc7[3,1] = 0;
dvdc7[4,1] = 0;
dvdc7[5,1] = 0;
dvdc7[6,1] = 0;
dvdc7[7,1] = 0;
dvdc7[8,1] = 0;
dvdc7[9,1] = -1;
dvdc7[10,1] = 0;
dvdc7[11,1] = 0;
dvdc7[12,1] = 0;
dvdc7[13,1] = -(a[1,1]^2)*a[2,1] -a[3,1];
dvdc7[14,1] = -4*a[1,1]*a[2,1];
dvdc7[15,1] = 0;
dvdc7[16,1] = 0;
dvdc7[17,1] = 0;
dvdc7[18,1] = 0;
dvdc7[19,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc7[20,1] = -4*(a[1,1]^2)*a[5,1];
dvdc7[21,1] = -10*a[1,1]*a[5,1];

dvdc7[22,1] = 0;
dvdc7[23,1] = -10*a[2,1];
dvdc7[24,1] = 0;
dvdc7[25,1] = 0;
dvdc7[26,1] = -(a[1,1]^4)*a[6,1] - 6*(a[1,1]^2)*a[2,1]*a[3,1] - a[10,1];
dvdc7[27,1] = -4*(a[1,1]^3)*a[6,1] - 12*a[1,1]*a[2,1]*a[3,1];
dvdc7[28,1] = -10*(a[1,1]^2)*a[6,1] - 10*a[3,1]*a[2,1];
dvdc7[29,1] = -20*a[1,1]*a[6,1];



dvdc8[1,1] = 0;
dvdc8[2,1] = 0;
dvdc8[3,1] = 0;
dvdc8[4,1] = 0;
dvdc8[5,1] = 0;
dvdc8[6,1] = 0;
dvdc8[7,1] = 0;
dvdc8[8,1] = 0;
dvdc8[9,1] = 0;
dvdc8[10,1] = -1;
dvdc8[11,1] = -4*a[1,1]*a[2,1];
dvdc8[12,1] = -( a[2,1] + a[4,1] );
dvdc8[13,1] = 0;
dvdc8[14,1] = 0;
dvdc8[15,1] = 0;
dvdc8[16,1] = 0;
dvdc8[17,1] = -10*(a[1,1]^2)*a[5,1];
dvdc8[18,1] = -4*a[1,1]*a[5,1];
dvdc8[19,1] = -( a[5,1] + a[7,1] );
dvdc8[20,1] = 0;
dvdc8[21,1] = 0;

dvdc8[22,1] = -10*(a[1,1]^2)*a[2,1];
dvdc8[23,1] = 0;
dvdc8[24,1] = -20*(a[1,1]^3)*a[6,1];
dvdc8[25,1] = -10*(a[1,1]^2)*(a[6,1] + a[2,1]*a[4,1]);
dvdc8[26,1] = -4*a[1,1]*(a[6,1] + 3*a[2,1]*a[4,1]);
dvdc8[27,1] = -(a[6,1] + 6*a[2,1]*a[4,1] + a[11,1]);
dvdc8[28,1] = 0;
dvdc8[29,1] = 0;



dvdc9[1,1] = 0;
dvdc9[2,1] = 0;
dvdc9[3,1] = 0;
dvdc9[4,1] = 0;
dvdc9[5,1] = 0;
dvdc9[6,1] = 0;
dvdc9[7,1] = 0;
dvdc9[8,1] = 0;
dvdc9[9,1] = 0;
dvdc9[10,1] = 0;
dvdc9[11,1] = -(a[1,1]^4);
dvdc9[12,1] = -(a[1,1]^3);
dvdc9[13,1] = -(a[1,1]^2);
dvdc9[14,1] = -a[1,1];
dvdc9[15,1] = 0;
dvdc9[16,1] = 0;
dvdc9[17,1] = 0;
dvdc9[18,1] = 0;
dvdc9[19,1] = 0;
dvdc9[20,1] = 0;
dvdc9[21,1] = 0;

dvdc9[22,1] = -(a[1,1]^5);
dvdc9[23,1] = -1;
dvdc9[24,1] = -15*(a[1,1]^4)*a[3,1];
dvdc9[25,1] = -(a[1,1]^5)*a[4,1] - 10*(a[1,1]^3)*a[3,1];
dvdc9[26,1] = -3*(a[1,1]^4)*a[4,1] - 6*(a[1,1]^2)*a[3,1];
dvdc9[27,1] = -6*(a[1,1]^3)*a[4,1] - 3*a[1,1]*a[3,1];
dvdc9[28,1] = -10*(a[1,1]^2)*a[4,1] - a[3,1];
dvdc9[29,1] = -15*a[1,1]*a[4,1];



dvdc10[1,1] = 0;
dvdc10[2,1] = 0;
dvdc10[3,1] = 0;
dvdc10[4,1] = 0;
dvdc10[5,1] = 0;
dvdc10[6,1] = 0;
dvdc10[7,1] = 0;
dvdc10[8,1] = 0;
dvdc10[9,1] = 0;
dvdc10[10,1] = 0;
dvdc10[11,1] = 0;
dvdc10[12,1] = 0;
dvdc10[13,1] = 0;
dvdc10[14,1] = 0;
dvdc10[15,1] = -1;
dvdc10[16,1] = 0;
dvdc10[17,1] = -5*a[1,1]*a[2,1];
dvdc10[18,1] = -( a[2,1] + a[4,1] );
dvdc10[19,1] = 0;
dvdc10[20,1] = 0;
dvdc10[21,1] = 0;

dvdc10[22,1] = 0;
dvdc10[23,1] = 0;
dvdc10[24,1] = -15*(a[1,1]^2)*a[5,1];
dvdc10[25,1] = -5*a[1,1]*a[5,1];
dvdc10[26,1] = -(a[5,1] + a[7,1]);
dvdc10[27,1] = 0;
dvdc10[28,1] = 0;
dvdc10[29,1] = 0;



dvdc11[1,1] = 0;
dvdc11[2,1] = 0;
dvdc11[3,1] = 0;
dvdc11[4,1] = 0;
dvdc11[5,1] = 0;
dvdc11[6,1] = 0;
dvdc11[7,1] = 0;
dvdc11[8,1] = 0;
dvdc11[9,1] = 0;
dvdc11[10,1] = 0;
dvdc11[11,1] = 0;
dvdc11[12,1] = 0;
dvdc11[13,1] = 0;
dvdc11[14,1] = 0;
dvdc11[15,1] = 0;
dvdc11[16,1] = -1;
dvdc11[17,1] = 0;
dvdc11[18,1] = 0;
dvdc11[19,1] = 0;
dvdc11[20,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc11[21,1] = -5*a[1,1]*a[2,1];

dvdc11[22,1] = 0;
dvdc11[23,1] = 0;
dvdc11[24,1] = 0;
dvdc11[25,1] = 0;
dvdc11[26,1] = 0;
dvdc11[27,1] = -(a[1,1]^3)*a[5,1] - a[8,1];
dvdc11[28,1] = -5*(a[1,1]^2)*a[5,1];
dvdc11[29,1] = -15*a[1,1]*a[5,1];



dvdc12[1,1] = 0;
dvdc12[2,1] = 0;
dvdc12[3,1] = 0;
dvdc12[4,1] = 0;
dvdc12[5,1] = 0;
dvdc12[6,1] = 0;
dvdc12[7,1] = 0;
dvdc12[8,1] = 0;
dvdc12[9,1] = 0;
dvdc12[10,1] = 0;
dvdc12[11,1] = 0;
dvdc12[12,1] = 0;
dvdc12[13,1] = 0;
dvdc12[14,1] = 0;
dvdc12[15,1] = 0;
dvdc12[16,1] = 0;
dvdc12[17,1] = -(a[1,1]^5);
dvdc12[18,1] = -(a[1,1]^4);
dvdc12[19,1] = -(a[1,1]^3);
dvdc12[20,1] = -(a[1,1]^2);
dvdc12[21,1] = -a[1,1];

dvdc12[22,1] = 0;
dvdc12[23,1] = 0;
dvdc12[24,1] = 0;
dvdc12[25,1] = 0;
dvdc12[26,1] = 0;
dvdc12[27,1] = 0;
dvdc12[28,1] = 0;
dvdc12[29,1] = 0;



dvdc13[1,1] = 0;
dvdc13[2,1] = 0;
dvdc13[3,1] = 0;
dvdc13[4,1] = 0;
dvdc13[5,1] = 0;
dvdc13[6,1] = 0;
dvdc13[7,1] = 0;
dvdc13[8,1] = 0;
dvdc13[9,1] = 0;
dvdc13[10,1] = 0;
dvdc13[11,1] = 0;
dvdc13[12,1] = 0;
dvdc13[13,1] = 0;
dvdc13[14,1] = 0;
dvdc13[15,1] = 0;
dvdc13[16,1] = 0;
dvdc13[17,1] = 0;
dvdc13[18,1] = 0;
dvdc13[19,1] = 0;
dvdc13[20,1] = 0;
dvdc13[21,1] = 0;
dvdc13[22,1] = -1;
dvdc13[23,1] = 0;
dvdc13[24,1] = -6*a[1,1]*a[2,1];
dvdc13[25,1] = -(a[2,1] + a[4,1]);
dvdc13[26,1] = 0;
dvdc13[27,1] = 0;
dvdc13[28,1] = 0;
dvdc13[29,1] = 0;



dvdc14[1,1] = 0;
dvdc14[2,1] = 0;
dvdc14[3,1] = 0;
dvdc14[4,1] = 0;
dvdc14[5,1] = 0;
dvdc14[6,1] = 0;
dvdc14[7,1] = 0;
dvdc14[8,1] = 0;
dvdc14[9,1] = 0;
dvdc14[10,1] = 0;
dvdc14[11,1] = 0;
dvdc14[12,1] = 0;
dvdc14[13,1] = 0;
dvdc14[14,1] = 0;
dvdc14[15,1] = 0;
dvdc14[16,1] = 0;
dvdc14[17,1] = 0;
dvdc14[18,1] = 0;
dvdc14[19,1] = 0;
dvdc14[20,1] = 0;
dvdc14[21,1] = 0;
dvdc14[22,1] = 0;
dvdc14[23,1] = -1;
dvdc14[24,1] = 0;
dvdc14[25,1] = 0;
dvdc14[26,1] = 0;
dvdc14[27,1] = 0;
dvdc14[28,1] = -(a[1,1]^2)*a[2,1] - a[3,1];
dvdc14[29,1] = -6*a[1,1]*a[2,1];



dvdc15[1,1] = 0;
dvdc15[2,1] = 0;
dvdc15[3,1] = 0;
dvdc15[4,1] = 0;
dvdc15[5,1] = 0;
dvdc15[6,1] = 0;
dvdc15[7,1] = 0;
dvdc15[8,1] = 0;
dvdc15[9,1] = 0;
dvdc15[10,1] = 0;
dvdc15[11,1] = 0;
dvdc15[12,1] = 0;
dvdc15[13,1] = 0;
dvdc15[14,1] = 0;
dvdc15[15,1] = 0;
dvdc15[16,1] = 0;
dvdc15[17,1] = 0;
dvdc15[18,1] = 0;
dvdc15[19,1] = 0;
dvdc15[20,1] = 0;
dvdc15[21,1] = 0;
dvdc15[22,1] = 0;
dvdc15[23,1] = 0;
dvdc15[24,1] = -(a[1,1]^6);
dvdc15[25,1] = -(a[1,1]^5);
dvdc15[26,1] = -(a[1,1]^4);
dvdc15[27,1] = -(a[1,1]^3);
dvdc15[28,1] = -(a[1,1]^2);
dvdc15[29,1] = -a[1,1];

retp(dvdc1~dvdc2~dvdc3~dvdc4~dvdc5~dvdc6~dvdc7~dvdc8~dvdc9
          ~dvdc10~dvdc11~dvdc12~dvdc13~dvdc14~dvdc15);

clear dvdc1, dvdc2, dvdc3, dvdc4, dvdc5, dvdc6, dvdc7, dvdc8, dvdc9,
      dvdc10, dvdc11, dvdc12, dvdc13, dvdc14, dvdc15;

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
      s_f1 = deff(s_c1,Emom,0);
    lc1 = s_f1'w*s_f1;
    clear s_f1;
  do until s_itr > maxsqez;

    s_c2=s_c-s_lm*s_dc;
      s_f2 = deff(s_c2,Emom,0);
    lc2 = s_f2'w*s_f2;
    clear s_f2;

    if lc1 <= lc2 and lc1 <= obj;

       retp(s_c1,s_itr-1); goto eoproc;

    else;

       s_c1=s_c2; s_lm=s_lm/2; lc1=lc2;
       s_itr=s_itr+1;

    endif;

  endo;
retp(s_c2,s_itr-1);

eoproc:
endp;


@  ------------------------------------------------------- @
@ Subroutine to compute squeezes.

  This subroutine compares the values of the objective function at the
  points c+s*dc and c+0.5*s*dc, with dc = the proposed change in the vector
  of parameters, and with step length s initially at 1. s is halved until
  minus the objective function stops declining.
@
proc (2) = sq(s_c,s_dc);
local s_c1,s_lm,s_itr,lc1,s_c2,lc2,s_f1,s_f2;


    s_c1=s_c - s_dc; s_lm=1/2; s_itr=1;
      s_f1 = saveit[1,idx]' - s_c1;
    lc1 = s_f1'w*s_f1;
    clear s_f1;
  do until s_itr > maxsqez;

    s_c2=s_c-s_lm*s_dc;
      s_f2 = saveit[1,idx]' - s_c2;
    lc2 = s_f2'w*s_f2;
    clear s_f2;

    if lc1 <= lc2 and lc1 <= obj;

       retp(s_c1,s_itr-1); goto eproc;

    else;

       s_c1=s_c2; s_lm=s_lm/2; lc1=lc2;
       s_itr=s_itr+1;

    endif;

  endo;
retp(s_c2,s_itr-1);

eproc:
endp;

@----------Subroutine to trap inversion errors-----------------@

proc invit(w);
local onemore,ww;
onemore = 1;
winv:
trap 1;
  ww = invpd(w);
trap 0;
  if scalerr(ww) > 0;
      if onemore == 1;
        call sysstate(14,1e-128);
        onemore = 0;
        goto winv;
      else;
        ww = -999999;
      endif;
  endif;

retp(ww);

endp;

@----------Subroutine to trap solving errors-----------------@

proc solvit(a,b);
local dc,onemore;
    onemore = 1;
    solveit:
    trap 1;
    dc = solpd(a,b);
    trap 0;
    if scalerr(dc) > 0;
      if onemore == 1;
        call sysstate(14,1e-128);
        onemore = 0;
        goto solveit;
      else;
        dc = -999999;
      endif;
    endif;

retp(dc);
endp;

@------------------------------------------------------------------------@

eop:

"Start time:  ";; timestr(ttt);
"End time:    ";; timestr(0);
"Start date:  ";; datestr(ddd);
"End date:    ";; datestr(0);

output off;
end;
