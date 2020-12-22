* Jonathan Norris
* Version 1.0.0 
*
* wildrw implements the wild cluster bootstrap with the Ramona Wolf 
* adjustment for multiple hypothesis testing 

* currently works with areg and nbreg. to use nbreg, must specify the method option.
* method() must list areg or nbreg in order of the outcome it goes with
* following the list of outcomes in the syntax

*********************************************************
* Steps *
* get sim t-values (tsim_s) on x for each outcome

* sort real t by largest t and start at s=1
* define maxt vector where each row takes the max across each row of (tsim_1 ... tsim_S) foreach s
* padj1 = #(maxt >= real t) / (reps)

* go to s=2 (2nd largest real t)
* define maxt vector where each row takes the max across each row of (tsim_2 ... tsim_S)
* e.g. tsim1 is dropped
* pint =  #(maxt >= real t) / (reps)
* padj2 = max(padj1, pint)

* repeat through s=S
*********************************************************

capture program drop wildrw
{
	program define wildrw, rclass
	syntax varlist [if] [in], ADJvars(varlist) /* 
		*/ MCONtrols1(varlist fv) /* 
		*/ Absorb(varlist) /* 
		*/ CLUster(varlist) /*
		*/ REPS(numlist) /*
		*/ [mcontrols2(varlist fv) mcontrols3(varlist fv) mcontrols4(varlist fv) /*
		*/ seed(numlist) /*
		*/ method(namelist) /*
		*/ sampby(varlist max=1)] /*
		*/
	
	preserve 
	marksample touse
	set seed `seed'
	local numYvars: word count `varlist'
	local numY = `numYvars'
	
	tokenize `adjvars'
	local adj "`1'"

	tokenize `absorb'
	local fe "`*'"
	tokenize `cluster'
	local clu "`*'"

	* setup local to denote what method goes with what var
	* in syntax method list must be in same order as outcome list
	local vpos = 1
	foreach m of local method {
		if "`method'" == "" {
			local method_`vpos' = "`vpos'_areg"
		}
		if "`m'" == "areg" {
			local method_`vpos' = "`vpos'_areg"
		}
		if "`m'" == "nbreg" {
			local method_`vpos' = "`vpos'_nbreg"   // only with nbreg in mind at this time
		}
		local ++vpos 
	}
	

	* count number of specs per outcome 
	local numSpec
	forvalues mcon = 1/4 {
		if "`mcontrols`mcon''" != "" {
			local numSpec = `numSpec' + 1
		}
	}

	* Wild Bootstrap *
	tempname b se 
	if "`sampby'" == "" {
		mat `b' = J(`numY'*`numSpec', 1, .)
		mat `se' = J(`numY'*`numSpec',1,.)
	}
	else {
		// x2 b/c running regs overs each group (e.g. female/male)
		mat `b' = J(`numY'*`numSpec'*2, 1, .)
		mat `se' = J(`numY'*`numSpec'*2,1,.)
	}
	local i = 1
	local S = 0
	local vpos = 1
	* initial varlist and tsim varlist
	local trealvar1
	local tsimvar1
	local jorder1

	foreach var of varlist `varlist' {
		foreach spec of numlist 1/`numSpec' {
			if "`sampby'" == "" {
				// run either areg or nbreg
				if "`method_`vpos''" == "`vpos'_areg" {
					// base regression
					qui areg `var' `adj'      ///
						`mcontrols`spec''  ///
						if `var' !=., absorb(`fe') cluster(`clu')
					mat `b'[`i',1] = _b[`adj']
					mat `se'[`i',1] = _se[`adj']
				}
				if "`method_`vpos''" == "`vpos'_nbreg" {
					qui nbreg `var' `adj'    ///
						`mcontrols`spec''    ///
						i.`fe' if `var' != ., vce(cluster `clu')
					mat `b'[`i',1] = _b[`adj']
					mat `se'[`i',1] = _se[`adj']
				}
			}
			else {
				// run either areg or nbreg
				if "`method_`vpos''" == "`vpos'_areg" {
					// base regression for sampby==1
					qui areg `var' `adj'      ///
						`mcontrols`spec''  ///
						if `var' !=. & `sampby' == 1,    ///
						absorb(`fe') cluster(`clu')
					mat `b'[`i',1] = _b[`adj']
					mat `se'[`i',1] = _se[`adj']

					// store actual t-value 
					tempvar tval_`var'_1_m`spec'
					cap drop `tval_`var'_1_m`spec''
					qui gen `tval_`var'_1_m`spec'' =  abs(`b'[`i',1] / `se'[`i',1])
					// store in local for cal in next step 
					local trealvar1 `trealvar1' `tval_`var'_1_m`spec''
			
					// wild cluster bootstrap
					qui boottest `adj', reps(`reps') svmat nograph seed(`seed')
					//mat `wildp'[`i',1] = r(p)    // wild cluster p-value 
					local wildp_`var'_1_m`spec' = r(p)    // wild cluster p-value 
					tempname tsim_`var'_1_m`spec'
					mat `tsim_`var'_1_m`spec'' = r(dist)
			
					// get abs of tsim vector & send to stata as a variable 
					tempname tsim
					mata: `tsim' = st_matrix("`tsim_`var'_1_m`spec''");
					mata: `tsim' = abs(`tsim');
					tempvar tsim_`var'_1_m`spec'
					getmata `tsim_`var'_1_m`spec''=`tsim', replace force
					local tsimvar1 `tsimvar1' `tsim_`var'_1_m`spec''
					
					local jorder1 `jorder1' `var'_1_m`spec'
		
					//local ++spec 
					local ++i 
					local ++S  // counter for how many hypothses to adjust

					// base regression for sampby==0
					qui areg `var' `adj'      ///
						`mcontrols`spec''  ///
						if `var' !=. & `sampby' == 0,    ///
						absorb(`fe') cluster(`clu')
					mat `b'[`i',1] = _b[`adj']
					mat `se'[`i',1] = _se[`adj']

					// store actual t-value 
					tempvar tval_`var'_0_m`spec'
					cap drop `tval_`var'_0_m`spec''
					qui gen `tval_`var'_0_m`spec'' =  abs(`b'[`i',1] / `se'[`i',1])
					// store in local for cal in next step 
					local trealvar1 `trealvar1' `tval_`var'_0_m`spec''
			
					// wild cluster bootstrap
					qui boottest `adj', reps(`reps') svmat nograph seed(`seed')
					//mat `wildp'[`i',1] = r(p)    // wild cluster p-value 
					local wildp_`var'_0_m`spec' = r(p)    // wild cluster p-value 
					tempname tsim_`var'_0_m`spec'
					mat `tsim_`var'_0_m`spec'' = r(dist)
			
					// get abs of tsim vector & send to stata as a variable 
					tempname tsim
					mata: `tsim' = st_matrix("`tsim_`var'_0_m`spec''");
					mata: `tsim' = abs(`tsim');
					tempvar tsim_`var'_0_m`spec'
					getmata `tsim_`var'_0_m`spec''=`tsim', replace force
					local tsimvar1 `tsimvar1' `tsim_`var'_0_m`spec''
					
					local jorder1 `jorder1' `var'_0_m`spec'
		
					//local ++spec 
					local ++i 
					local ++S  // counter for how many hypothses to adjust
				}
				if "`method_`vpos''" == "`vpos'_nbreg" {
					// for sampby==1
					qui nbreg `var' `adj'    ///
						`mcontrols`spec'' i.`fe'   /// 
						if `var' != . & `sampby' == 1,   ///
						vce(cluster `clu')
					mat `b'[`i',1] = _b[`adj']
					mat `se'[`i',1] = _se[`adj']

					// store actual t-value 
					tempvar tval_`var'_1_m`spec'
					cap drop `tval_`var'_1_m`spec''
					qui gen `tval_`var'_1_m`spec'' =  abs(`b'[`i',1] / `se'[`i',1])
					// store in local for cal in next step 
					local trealvar1 `trealvar1' `tval_`var'_1_m`spec''
			
					// wild cluster bootstrap
					qui boottest `adj', reps(`reps') svmat nograph seed(`seed')
					//mat `wildp'[`i',1] = r(p)    // wild cluster p-value 
					local wildp_`var'_1_m`spec' = r(p)    // wild cluster p-value 
					tempname tsim_`var'_1_m`spec'
					mat `tsim_`var'_1_m`spec'' = r(dist)
			
					// get abs of tsim vector & send to stata as a variable 
					tempname tsim
					mata: `tsim' = st_matrix("`tsim_`var'_1_m`spec''");
					mata: `tsim' = abs(`tsim');
					tempvar tsim_`var'_1_m`spec'
					getmata `tsim_`var'_1_m`spec''=`tsim', replace force
					local tsimvar1 `tsimvar1' `tsim_`var'_1_m`spec''
					
					local jorder1 `jorder1' `var'_1_m`spec'
		
					//local ++spec 
					local ++i 
					local ++S  // counter for how many hypothses to adjust

				// for sampby==0
					qui nbreg `var' `adj'    ///
						`mcontrols`spec'' i.`fe'   /// 
						if `var' != . & `sampby' == 0,   ///
						vce(cluster `clu')
					mat `b'[`i',1] = _b[`adj']
					mat `se'[`i',1] = _se[`adj']

					// store actual t-value 
					tempvar tval_`var'_0_m`spec'
					cap drop `tval_`var'_0_m`spec''
					qui gen `tval_`var'_0_m`spec'' =  abs(`b'[`i',1] / `se'[`i',1])
					// store in local for cal in next step 
					local trealvar1 `trealvar1' `tval_`var'_0_m`spec''
			
					// wild cluster bootstrap
					qui boottest `adj', reps(`reps') svmat nograph seed(`seed')
					//mat `wildp'[`i',1] = r(p)    // wild cluster p-value 
					local wildp_`var'_0_m`spec' = r(p)    // wild cluster p-value 
					tempname tsim_`var'_0_m`spec'
					mat `tsim_`var'_0_m`spec'' = r(dist)
			
					// get abs of tsim vector & send to stata as a variable 
					tempname tsim
					mata: `tsim' = st_matrix("`tsim_`var'_0_m`spec''");
					mata: `tsim' = abs(`tsim');
					tempvar tsim_`var'_0_m`spec'
					getmata `tsim_`var'_0_m`spec''=`tsim', replace force
					local tsimvar1 `tsimvar1' `tsim_`var'_0_m`spec''
					
					local jorder1 `jorder1' `var'_0_m`spec'
		
					//local ++spec 
					local ++i 
					local ++S  // counter for how many hypothses to adjust
				}
			}	

			if "`sampby'" == "" {
				// store actual t-value 
				tempvar tval_`var'_m`spec'
				cap drop `tval_`var'_m`spec''
				qui gen `tval_`var'_m`spec'' =  abs(`b'[`i',1] / `se'[`i',1])
				// store in local for cal in next step 
				local trealvar1 `trealvar1' `tval_`var'_m`spec''
	
				// wild cluster bootstrap
				qui boottest `adj', reps(`reps') svmat nograph seed(`seed')
				//mat `wildp'[`i',1] = r(p)    // wild cluster p-value 
				local wildp_`var'_m`spec' = r(p)    // wild cluster p-value 
				tempname tsim_`var'_m`spec'
				mat `tsim_`var'_m`spec'' = r(dist)
		
				// get abs of tsim vector & send to stata as a variable 
				tempname tsim
				mata: `tsim' = st_matrix("`tsim_`var'_m`spec''");
				mata: `tsim' = abs(`tsim');
				tempvar tsim_`var'_m`spec'
				getmata `tsim_`var'_m`spec''=`tsim', replace force
				local tsimvar1 `tsimvar1' `tsim_`var'_m`spec''
				
				local jorder1 `jorder1' `var'_m`spec'
	
				//local ++spec 
				local ++i 
				local ++S  // counter for how many hypothses to adjust
			}
		}
		local ++vpos  // keeps track of position for method call
	}

	* order real t highest to lowest 
	* select largest sim t out of each rows sim t's iterating from j=1/S for j...S
	* note: recall `S' is defined in loop above & is number of hypotheses to adjust (num dep. vars)
	tempname wildp 
	if "`sampby'" == "" mat `wildp' = J(`numY'*`numSpec',1,.)
	if "`sampby'" != "" mat `wildp' = J(`numY'*`numSpec'*2,1,.)
	local varlist1 `varlist'
	local jorder   
	forvalues j = 1/`S' {
		tempvar realtmax`j'
		cap drop `realtmax`j''
		egen `realtmax`j'' = rowmax(`trealvar`j'')
		tempvar maxt`j'
		cap drop `maxt`j''
		egen `maxt`j'' = rowmax(`tsimvar`j'') 
		// initialize local names for next iteration
		local k = `j' + 1
		local jorder`k'
		local trealvar`k'
		local tsimvar`k'
		// setup varlists for next iteration 
		foreach v of local jorder`j' {
			if `tval_`v'' != `realtmax`j'' {
				local jorder`k' `jorder`k'' `v'
				local tsimvar`k' `tsimvar`k'' `tsim_`v''
				local trealvar`k' `trealvar`k'' `tval_`v''
			}
			if `tval_`v'' == `realtmax`j'' {
				tempvar tval_`j'_`v'
				cap drop `tval_`j'_`v''
				gen `tval_`j'_`v'' = `tval_`v''               // denote the j order of the true t 
				local jorder `jorder' `v'                     // put varnames in order of highest to lowest real t

				mat `wildp'[`j',1] =  `wildp_`v''
			}
		}
	}

	* get RW adjusted p-values for j=1 to S
	//forvalues j = 1/`S' {
		//foreach spec of numlist 1/`numSpec' {
	local j = 1
	foreach v of local jorder {
		if `j' == 1 {
			count if `maxt`j'' >= `tval_`j'_`v'' & `maxt`j'' != .
			local cnum1 = r(N)
			tempvar prw_`j'_`v'
			cap drop `prw_`j'_`v''
			gen `prw_`j'_`v'' = `cnum1' / `reps'
		}
		else {
			qui count if `maxt`j'' >= `tval_`j'_`v'' & `maxt`j'' != .
			local cnum = r(N)
			tempvar pint_`j'_`v'
			cap drop `pint_`j'_`v''
			gen `pint_`j'_`v'' = `cnum' / `reps'    // initial p-value 

			// for j=2 comparison below, first time need to define b/c stata can't see it
			tempvar prw_l_`v'
			gen `prw_l_`v'' = `cnum1' / `reps'  // j=1 prw 

			tempvar prw_`j'_`v' 
			egen `prw_`j'_`v'' = rowmax(`prw_l_`v'' `pint_`j'_`v'') if `prw_l_`v'' !=. & `pint_`j'_`v'' !=.   // enforce  monotonicity 

			// setup for next iteration 
			tempvar prw_l_`v'
			cap drop `prw_l_`v''
			gen `prw_l_`v'' = `cnum' / `reps'
		}
		local ++j
	}
	
	* now assign prw_`j'_m`spec' label for proper outcome it is attached to using the jorder local
	tempname prw 
	if "`sampby'" == "" mat `prw' = J(`numY'*`numSpec', 1, .)
	if "`sampby'" != "" mat `prw' = J(`numY'*`numSpec'*2, 1, .)
	local rowlabels
	local j = 1
	foreach v of local jorder {
		tempname p_`j'_`v'
		mkmat `prw_`j'_`v'' in 1, matrix(`p_`j'_`v'')
		mat `prw'[`j',1] = `p_`j'_`v''[1,1] 

		local rowlabels `rowlabels' `v'
		local ++j
	}

	mat rownames `prw' = `rowlabels'
	mat rownames `wildp' = `rowlabels'
	//mat list `prw' 

	foreach v of varlist `adj' {
		local j = 1
		dis "Independent variable: `adj'"
		dis "Outcome Variables: `rowlabels'"
		dis "Number of resamples: `reps'"
		dis _newline
		dis "{hline 78}"
		dis " Outcome Variable |  Wild Cluster p-value   Romano-Wolf p-value"
		dis "{hline 18}+{hline 57}"
		foreach var of local rowlabels {
			//foreach mcon of numlist 1/`numSpec' {
				local pvalwild = `wildp'[`j',1]
				local pvalrw = `prw'[`j',1]
            	display as text %17s abbrev("`var'",17) " {c |}     "   /*
            	*/ as result %6.4f `pvalwild'    "                " /*
            	*/ as result %6.4f `pvalrw'      "                "  

            	local ++j
            //}
		}
	}

	return matrix prw = `prw'
	return matrix wildp = `wildp'

	restore
end
}




