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
output file=ewdata.out reset;

@------Name of the Gauss data set.---------------------------------@

ttdata = "scomp";
ewdata = "ewdata";

@--------Open GAUSS data set. -------------------------------------@

    closeall f1,f2;
    open f1=^ttdata;

@------------Create New GAUSS data set. -------------------------------@

let vnames = cusip sic year nk gk cf inv ik sk life brate divs sales;

create f2 = ^ewdata with ^vnames,0,8;

@----------OPEN ALL OF THE FILES---------------------------------------@

load streq;
  str = streq[.,2];
  equ = streq[.,3];
load taxes;
  t    = taxes[.,2];
  itcr = taxes[.,3];
load pki;
load ppi;

@-------------NEXT, CONSTRUCT THE ITC RATE----------------------------@

u = itcr.*str./(str + equ);

@--------------------REPEAT FOR EACH COMPANY-------------------------@

do until eof(f1);

 dta = readr(f1,1);                  $dta[1,1];
   nyr = dta[1,42];
   byr = dta[1,3];
   eyr = byr + nyr - 1;
 dta= dta|(readr(f1,nyr-1));

   invraw = dta[.,4];
   act    = dta[.,5];
   lct    = dta[.,6];
   gkraw  = dta[.,8];
   nkraw  = dta[.,9];
   sales  = dta[.,12];
   iex    = dta[.,15];
   opinc  = dta[.,13];
   depr   = dta[.,14];
   tax    = dta[.,16];
   inc    = dta[.,40];
   fgraw  = dta[.,39];
   i      = dta[.,35];
   sk     = missrv(dta[.,36],0);
   imeth  = missrv(dta[.,27],0);
   cdiv   = missrv(dta[.,18],99999);
   pdiv   = missrv(dta[.,17],99999);
   yr     = dta[.,3];   @ if dta[1,1]$=="00282410";  yr~depr~inc;end;   endif;   @
   cusip  = dta[1,1];
   icode  = dta[1,2];
   brate  = missrv(dta[.,41],9999);

@--------Screen for data on finished goods inventories.-------------@
   mm = miss(1,1);
   invok = fgraw .== mm;
   invok = 1 - (sumc(invok) > 0);

/*--------------------------------------------------------------------
*     The next section constructs two measures of the replacement
*     value of the capital stock.
*--------------------------------------------------------------------*/

    dind = depr .== 0;
    numiss = sumc(dind);
    ltt = gkraw[2:rows(gkraw),1] + i[1:(rows(gkraw)-1),1] - sk[1:(rows(gkraw)-1),1];
    ltt = ltt./(depr[2:rows(depr),1] + dind[2:rows(depr),1]).*(1-dind[2:rows(depr),1]);
    slt = sumc(ltt);
    l = slt/(nyr - 1 - numiss);

@--------------Next, define the capital stock iteratively-----------@
        nk = nkraw;
        gk = gkraw;
        nk[1,1] = nkraw[1,1];
        gk[1,1] = gkraw[1,1];
        isic = floor(dta[1,2]/100)-19;
        j=2;
        do while j <= nyr;
          nk[j,1] = (nk[j-1,1]*pki[j+byr-77,isic]/pki[j-1+byr-77,isic] + i[j,1] - sk[j,1])*(1 - 1/l);
          gk[j,1] = (gk[j-1,1]*pki[j+byr-77,isic]/pki[j-1+byr-77,isic] + i[j,1] - sk[j,1])*(1 - 1/l);
          j = j + 1;
        endo;
@-------------Divide by the relative price of capital goods----------@
          nk = nk./(pki[byr-77+1:eyr-77+1,isic]/ppi[byr-77+1:eyr-77+1,isic]);
          gk = gk./(pki[byr-77+1:eyr-77+1,isic]/ppi[byr-77+1:eyr-77+1,isic]);
          sk = sk./(pki[byr-77+1:eyr-77+1,isic]/ppi[byr-77+1:eyr-77+1,isic]);
          i  =  i./(pki[byr-77+1:eyr-77+1,isic]/ppi[byr-77+1:eyr-77+1,isic]);
/*
*--------------------------------------------------------------------
* this section constructs the replacement value of total inventories
*--------------------------------------------------------------------
*/

@----------------First, make the indicator for the inventory method.----@
         idx = imeth .> 2;
         imeth = imeth.*(1-idx) + idx;

@----------------first, the easy part--fifo--------------------------@
         invt = invraw;
@-------------next, the hard part--all or part lifo------------------@

@------------------all lifo------------------------------------------@
        if imeth == 2;
          jj = 2;
          do while jj <= nyr;
            if invraw[jj,1] >= invraw[jj-1,1];
              invt[jj,1] = invt[jj-1,1]*ppi[jj+byr-77,isic]/ppi[jj-1+byr-77,isic] +
                          invraw[jj,1] - invraw[jj-1,1];
            else;
              invt[jj,1] = (invt[jj-1,1] + invraw[jj,1] - invraw[jj-1,1])
                           *ppi[jj+byr-77,isic]/ppi[jj-1+byr-77,isic];
            endif;
          jj = jj + 1;
          endo;
        endif;
/*
*---------------------------------------------------------------------
*     this section constructs and spits the final regression variables
*---------------------------------------------------------------------
*/
        jj = 1;
        do while jj <= nyr;
          cf = depr[jj,1] + inc[jj,1];       yr[jj,1];;cf;
          divs = cdiv[jj,1] + pdiv[jj,1];
          z=writer(f2,(cusip~icode~yr[jj,1]~nk[jj,1]~gk[jj,1]~
                   cf~invt[jj,1]~i[jj,1]~sk[jj,1]~l~brate[jj,1]~divs~sales[jj,1]));
        jj = jj + 1;
        endo;
      endo;
      end;
