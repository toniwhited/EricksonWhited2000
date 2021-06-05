@ 1  cusip              @
@ 2  sic                @
@ 3  year               @
@ 4  invt               @
@ 5  act                @
@ 6  lct                @
@ 7  at                 @
@ 8  ppegt              @
@ 9  ppent              @
@10  dltt               @
@11  dltis              @
@12  sale               @
@13  oibdp              @
@14  dp                 @
@15  xint               @
@16  txt                @
@17  dvp                @
@18  dvc                @
@19  prcc               @
@20  csho               @
@21  dlc                @
@22  cogs               @
@23  xsga               @
@24  xlr                @
@25  xpr                @
@26  emp                @
@27  invval             @
@28  dd1                @
@29  dd2                @
@30  dd3                @
@31  dd4                @
@32  dd5                @
@33  sstk               @
@34  prstkc             @
@35  capx               @
@36  sppe               @
@37  aqc                @
@38  dltr               @
@39  fginv              @
@40  inc                @
@41  rate               @
@42  nyr                @

new;
output file=de.out reset;

@------Name of the Gauss data set.---------------------------------@

ttdata = "scomp";
dedata = "de";

@--------Open GAUSS data set. -------------------------------------@

    closeall f1,f2;
    open f1=^ttdata;

@------------Create New GAUSS data set. -------------------------------@

let vnames = cusip sic year tval eqty bookass bookdet;

create f2 = ^dedata with ^vnames,0,8;

@----------Open all of the files---------------------------------------@

      load interest;
      divyd = interest[.,3];
      baa   = interest[.,2];
      load matdist;
      matdist = matdist[.,2:cols(matdist)];

@-------Set the assumed maturity of long-term debt.--------------------@

k = 20;

@----------Construct the coupon and discount rates---------------------@

      disc = 1/(baa[21:rows(baa),1]/100 + 1);
      rr   = baa[21:rows(baa),1]/100;
      c = zeros(19,k);
      rate = zeros(19,k);
      jj = 1;
      do while jj <= k;
        c[.,jj] = (baa[(21+jj-k):(rows(baa)+jj-k),1]/100);
        rate[.,jj] = (baa[(21+jj-k):(rows(baa)+jj-k),1]./baa[21:rows(baa),1]);
        jj = jj + 1;
        endo;

@-------Read in the data.----------------------------------------------@

qq = 1;
do until eof(f1);


 dta = readr(f1,1);
 byr = dta[1,3];
 nyr = dta[1,42];
 eyr = byr + nyr - 1;
 dta= dta|(readr(f1,nyr-1));
 v = zeros(nyr,8);      $dta[1,1];


 v[.,1] = dta[.,10];
 v[.,2] = dta[.,15];
 v[.,3] = dta[.,21];
 v[.,4] = missrv(dta[.,28],0);
 v[.,5] = dta[.,29];
 v[.,6] = dta[.,30];
 v[.,7] = dta[.,31];
 v[.,8] = dta[.,32];

 w1 = dta[.,17];
 w2 = dta[.,19];
 w3 = dta[.,20];

@------Screen for data on d2-d5----------------------------------------@

  ms = miss(1,1);

  m1 = v[.,4] .== ms .or
       v[.,5] .== ms .or
       v[.,6] .== ms .or
       v[.,7] .== ms .or
       v[.,8] .== ms;


@-------------Calculate the market value of equity.---------------------@
          prstk = w1./(divyd[(dta[1,3]-56):(dta[1,3]-56+dta[1,42]-1)]./100);
          cmstk = w2.*w3;
          eqty  = prstk + cmstk;

@-------This section constructs the initial bsw series--------------@

  d=zeros(nyr,k);  @Columns are years, rows are maturity. @
  dish = zeros(nyr,1);

/* Set the maturity distribution in the first year equal to a uniform. */
  matdist = 0.05*ones(1,20);
  d[1,.] = (matdist)*(v[1,1]+v[1,4]);

/* Set debt issuance in the first year equal to the amount assumed due in k years */

  dish[1,1] = d[1,k];

/* Now update the rest of the years using information on new debt issuance. */

jj = 2;
do while jj <= nyr;

    del = v[jj,1] + v[jj,4] - (v[(jj-1),1] + v[(jj-1),4]);
    dish[jj,1] = del;
    ret = d[(jj-1),1];
    d[jj,1:(k-1)] = d[jj-1,2:k];
    if (del + ret) > 0;
      d[jj,k] = del + ret;
    else;
      d[jj,k] = 0;
      denom = v[(jj-1),1] + v[(jj-1),4] - ret;
      if denom /= 0;
        pcg = (v[jj,1]+v[jj,4])/(v[(jj-1),1]+v[(jj-1),4] - ret);
        d[jj,1:(k-1)] = d[jj,1:(k-1)]*pcg;
      endif;
    endif;

jj = jj + 1;
endo;

/*----------------------------------------------------------------------
*   This section modifies the bsw calculations to make use of the data
*   on debt due in one to five years
------------------------------------------------------------------------*/

jj = 1;
do while jj <= nyr;

@------If the data exist, replace the calculated d's with the actual.---@
  if m1[jj,1] == 0;
      diff = d[jj,1:5] - v[jj,4:8];
      tdiff = sumc(diff');
      d[jj,1:5] = v[jj,4:8];
@---------Construct very long term debt--------------------------------@
      tvltd = sumc(d[jj,6:k]');
@--------Initialize the "trial" debt values----------------------------@
      dd=zeros(1,k-5);
@----Distribute the difference proportionally over the remaining years-@
      if tvltd > 0 ;
        pcg = d[jj,6:k]/tvltd;
        dd = d[jj,6:k] + pcg*tdiff;
@-----If dd(j) is negative, set it to zero and update the totals-------@
        isneg = (dd .< 0);
        sumneg = abs(sumc((isneg.*dd)'));
        if (cols(isneg) - sumc(isneg')) == 0;
        dd = zeros(1,k-5);
        else;
        dd = dd.*(1 - isneg)*sumneg/(cols(isneg) - sumc(isneg'));
        endif;
       endif;
        d[jj,6:k] = dd;
   endif;

jj = jj + 1;
endo;

/*-----------------------------------------------------------------------
*  This section adjusts the above series for interest expense and uses
*  it to calculate market value
*------------------------------------------------------------------------*/

note = zeros(nyr,1);
jj = 1;
do while jj <= nyr;

  note[jj,1] = v[jj,3] - v[jj,4];
  tltd = sumc(d[jj,.]');
  alldet = tltd + note[jj,1];
  if alldet /= 0;
@--------Calculate ratios of debt due in mm years to total debt--------@
    r = d[jj,.]/alldet;
    aanote = note[jj,1]/(tltd + note[jj,1]);
@----Calculate interest expense divided by total debt------------------@
    sumrb = aanote*rr[byr-77+jj,1] + sumc((r.*c[byr-77+jj,.])');
@-------------------"Calculate" total debt-----------------------------@
    tbookd = v[jj,2]/sumrb;
@-----------------"Recalculate" the d's--------------------------------@
    d[jj,.] = r*tbookd;
    note[jj,1] = aanote*tbookd;
  endif;

jj = jj + 1;
endo;


@---------Do the market value calculations-----------------------------@

jj = 1;
do while jj <= nyr;
  tval = note[jj,1] + sumc( ((disc[jj,1]^seqa(1,1,k)).*(1-rate[jj,.]') + rate[jj,.]').*d[jj,.]');
  alldet = v[jj,1] + v[jj,3] - v[jj,4];
  if tval <0; qq; tval = v[jj,1]; endif;
  zz=writer(f2,(dta[jj,1:3]~tval~eqty[jj,1]~dta[jj,7]~alldet));
jj = jj + 1;
endo;
qq = qq + 1;
endo;
end;
