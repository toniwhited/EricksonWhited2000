new;
ttt = time;
ddd = date;

@------Name of the Gauss data set.---------------------------------@

indata    = "ew2000rawdata";
outdata   = "scomp";

output file= screen.out reset;

@--------Open GAUSS data set. -------------------------------------@

    closeall f1,f2;
    open f1=^indata;

@ 1 cusip    @
@ 2 sic      @
@ 3 year     @
@ 4 invt     @
@ 5 act      @
@ 6 lct      @
@ 7 at       @
@ 8 ppegt    @
@ 9 ppent    @
@10 dltt     @
@11 dltis    @
@12 sale     @
@13 oibdp    @
@14 dp       @
@15 xint     @
@16 txt      @
@17 dvp      @
@18 dvc      @
@19 prcc     @
@20 csho     @
@21 dlc      @
@22 cogs     @
@23 xsga     @
@24 xlr      @
@25 xpr      @
@26 emp      @
@27 invval   @
@28 dd1      @
@29 dd2      @
@30 dd3      @
@31 dd4      @
@32 dd5      @
@33 sstk     @
@34 prstkc   @
@35 capx     @
@36 sppe     @
@37 aqc      @
@38 dltr     @
@39 fginv    @
@40 inc      @
@41 rate     @
@42 nyr      @


vnames = getname(indata)|"nyr";


create f2 = ^outdata with ^vnames,0,8;

nyr = 19;
nco = rowsf(f1)/nyr;
yr = 5;

nfirm = 0;
i = 1;
do while i <=nco;

output off;
locate 3,1;
if (i/100)==floor(i/100);
"Firm Number " i;
endif;
output on;

x=readr(f1,nyr);

x[.,indcv("inc",vnames)] = rev(x[.,indcv("inc",vnames)]); /*This does the fix*/

/* Make indicators for the variables that you care about. */

mm = miss(1,1);
@
idx = (x[.,4  7 8 10 12 14 15 16 18 19 20 21  35 40] .== mm);
@
        let idx1 =   invt at ppegt dltt sale dp xint txt dvc prcc csho dlc capx inc;
idx = x[.,indcv(idx1,vnames)].==mm;
idx = sumc(idx');
idx = idx .== 0;


maxidx = 1;
jj = rows(idx)-yr;
do while jj >= 1;

if idx[jj,1] == 0;
maxidx = jj+1;
goto escape;
endif;

jj = jj - 1;
endo;
escape:

x=x[maxidx:rows(x),.];

/* This gets rid of firms with big mergers. */

bleh = 1 - (abs(missrv(x[.,37],0)./x[.,7]) .> 0.25);   @ 1 is good; 0 is bad. @

if sumc(bleh[(rows(x)-yr+1):rows(x),1]) == yr;

/* Find the lowest observation number such that everything is there. */

maxbleh = 1;
jj = rows(x)-yr+1;
do while jj >= 1;

if bleh[jj,1] == 0;
maxbleh = jj+1;
goto tonimw;
endif;

jj = jj - 1;
endo;
tonimw:

x=x[maxbleh:rows(x),.];
nyear = rows(x).*ones(rows(x),1);

x=x~nyear;

zz=writer(f2,x);
nfirm = nfirm + 1;
endif;

i = i + 1;
endo;

"Number of firms = " nfirm;

"Start time:  ";; timestr(ttt);
"End time:    ";; timestr(0);
"Start date:  ";; datestr(ddd);
"End date:    ";; datestr(0);

output off;
end;
