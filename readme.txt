To replicate the results, put all of the files in this zipped folder into the same directory.

Then run the programs in this order

Data manipulation programs:

1. screentt.g           This program gets rid of missing values and large acquisitions.
                        It also correctly sorts cash flow, firm by firm.
2. EWassets.g           This program makes the replacement values of the capital stock and inventories.
3. ewdebt.g             This program makes the market value of debt and equity.
4. dset.g               This program makes the final data set.

Estimation programs:

5. ewgmm.g              This program produces all of the results for tables 2-8.
6. idtest.g             This program produces all of the results for table 1.
7: muxmuy.g             This program produces all of the results for table 9.

Bootstrap programs:

8.  ewbootstrapBond.g        This program calculates the bootstrapped critical values for the bond interaction model.
9.  ewbootstrapSize.g        This program calculates the bootstrapped critical values for the size interaction model.
10. ewbootstrapBondSize.g    This program calculates the bootstrapped critical values for the bond and size interaction model.

---------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------
Data files:

EW2000rawdata.dat       A Gauss file containing the data pulled from the 1996 Compustat industrial files
                        for our original sample of firms.
divtax.fmt              Dividend tax rates
interest.fmt            Interest rates
matdist.fmt             Aggregate maturity distribution of debt
pki.fmt                 Index of capital goods prices
ppi.fmt                 Producer price index
streq.fmt               Mix of structures and equipment
taxes.fmt               Income tax rates
