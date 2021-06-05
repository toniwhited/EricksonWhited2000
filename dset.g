new;
output file=dset.out reset;

@--------------Number of observations------------------------------@

nyr = 4;

@------Name of the Gauss data set.---------------------------------@

ewdata   = "ewdata";
dedata   = "de";
dsetdata = "ewdset";

@--------Open GAUSS data set. -------------------------------------@

    closeall f1,f2,f3;
    open f1=^ewdata;
    open f2=^dedata;

@------------Create New GAUSS data set. -------------------------------@

let vnames = cusip sic year ik q cf k da rate payrat sales mtb qd qr assets;

create f3 = ^dsetdata with ^vnames,0,8;


@-------------Read in the little aggregate data files.-----------------@

      load interest;
      baa1   = interest[.,2];
      baa1 = baa1/100;

      load divtax;
      divtax1 = divtax[.,2];

      load taxes;
      cptax1 = taxes[.,2];

@-------------Read in the COMPUSTAT data.------------------------------@


y=readr(f1,rowsf(f1));
x=readr(f2,rowsf(f2));

@----------Set the number of cross-sectional observations.-----------@

nobs = sumc(y[.,3].==95-nyr+1);


@ewdata:    cusip sic year nk gk cf inv ik sk life brate divs sales@
@           1     2   3    4  5  6  7   8  9  10   11    12   13   @

@de:        cusip sic year tval eqty bookass bookdet@
@           1     2   3    4    5    6       7      @

@-------------Make the bitty data matrices the same as the big ones.---@


baa = zeros(rows(y),1);
divtax = zeros(rows(y),1);
cptax = zeros(rows(y),1);

for yy (77, 95, 1);

baa    = baa    + (y[.,3].==yy)*baa1[(yy-56),1];
divtax = divtax + (y[.,3].==yy)*divtax1[(yy-76),1];
cptax  = cptax  + (y[.,3].==yy)*cptax1[(yy-76),1];

endfor;

@----------------Create the variables.----------------------------------@

invt   = y[.,7];
tval   = x[.,4];
eqty   = x[.,5];
gk     = y[.,5];
at     = x[.,6];
bookd  = x[.,7];
life   = y[.,10];

qnum = (tval + eqty - invt);
qnum = qnum[1:(rows(qnum)-1),.];


z = (1./life)./((1./life) + (baa./0.95)) ;
xx = cptax.*z.*gk.*(1-divtax)./0.95 ;
qrnum = (1./(1 - cptax)).*(tval + eqty - xx - invt);
qdnum = (1./(1 - cptax)).*(tval + 0.95.*(eqty - xx)./(1-divtax) - invt);
qrnum = qrnum[1:(rows(qrnum)-1)];
qdnum = qdnum[1:(rows(qdnum)-1)];

gk   = gk[1:(rows(gk)-1)];

q    = qnum./gk;
qr   = qrnum./gk;
qd   = qdnum./gk;

bnum  = bookd + eqty;
@bnum  = eqty;@
bnum  = bnum[1:(rows(bnum)-1)];
at    = at[1:(rows(at)-1)];
bookd = bookd[1:(rows(bookd)-1)];
btm   = bnum./(gk+invt[1:(rows(invt)-1)]);


cusip  = y[2:rows(y),1];
sic    = y[2:rows(y),2];
year   = y[2:rows(y),3];
i      = y[2:rows(y),8];
sk     = y[2:rows(y),9];
rate   = y[2:rows(y),11];
cf     = y[2:rows(y),6];
divs   = y[2:rows(y),12];
sales  = y[2:rows(y),13];

i = (i - sk);
sk = sk;
ik = i./gk;
da = tval./(tval + eqty);
da = da[1:(rows(da)-1),1];

cf     = cf./gk;

sales  = sales./gk;

outmx = cusip~sic~year~ik~q~cf~gk~da~rate~divs~sales~btm~qd~qr~at;

ecrit = zeros(nyr*nobs,cols(outmx));



for qq (1, nyr, 1);

choix = (year .== qq+95-nyr);
choix = miss(choix,0);

choix = packr(outmx~choix);

choix = choix[.,1:(cols(choix)-1)];


ecrit[(qq-1)*nobs+1:qq*nobs,.] = choix;


endfor;



zz = writer(f3,ecrit);


end;
