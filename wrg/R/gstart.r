gstart<-function(x)  {		
## check basic format of x		
	if(class(x)!="data.frame") {stop("mle_gamma takes a structured dataframe input, use mleframe")}	
	if(ncol(x)!=3)  {stop("mle_gamma takes a structured dataframe input, use mleframe")}	
	xnames<-names(x)	
	if(xnames[1]!="left" || xnames[2]!="right"||xnames[3]!="qty")  {	
		 stop("mle_gamma takes a structured dataframe input, use mleframe") 
	}	
		
## disable warnings that likely will be issued for large data sets on the mlefit call due to Abpval		
	options(warn=-1)	
		
## subview the data according to failure or suspension		
	dafDF<- x[which(x$right>0),]	
	dasDF<- x[which(x$right<0),]	
		
# failure/suspension ratio		
	f2s<-sum(dafDF$qty)/ sum(dasDF$qty)	
## without consideration of late suspensions (of which there are none in the challenge case)		
	if(f2s >.95)  {	mat_line<-1
	}else if(f2s > .85)  {	mat_line<-2
	}else if(f2s > .75)  {	mat_line<-3
	}else if(f2s > .65)  {	mat_line<-4
	}else if(f2s > .55)  {	mat_line<-5
	}else if(f2s > .45)  {	mat_line<-6
	}else if(f2s > .35)  {	mat_line<-7
	}else if(f2s > .25)  {	mat_line<-8
	}else if(f2s > .15)  {	mat_line<-9
	}else if(f2s > .075)  {	mat_line<-10
	}else  	mat_line<-11
		
## usng the sub view of  failures (as defined by intervals)		
## to get the  daf_mean		
	row.sums<-apply(dafDF[,1:2], 1,sum)	
	calcDF<-data.frame(row.sums, dafDF$qty, .5)	
	row.prods<-apply(calcDF,1,prod)	
	fail_data_mean<-sum(row.prods)/sum(dafDF$qty)	
		
## get the weibull beta as a key to weibull 2 gamma correlation		
	wfit<-mlefit(x)	
	beta.w<-unname(wfit[2])	
		
M1a	<- matrix(	
c(0.04079731, 0.3058897, 0.6589425,		
0.04674198, 0.3109975, 0.6447629,		
0.04266909, 0.3583335, 0.6042592,		
0.03730037, 0.4039402, 0.5618611,		
0.04218690, 0.3991795, 0.5672162,		
0.03511832, 0.4537277, 0.5177023,		
0.03810429, 0.4626558, 0.5031752,		
0.03323306, 0.5176418, 0.4512175,		
0.03606437, 0.4872300, 0.4809901,		
0.03973339, 0.4976159, 0.4697456,		
0.04117619, 0.5010652, 0.4604034),		
nrow=10, ncol=3)		
		
M1b	<- matrix(	
c(0.9627292, -0.7484642, 1.0320645,		
0.3381235, -0.3209939, 0.9335946,		
0.5914780, -0.4527554, 0.9130913,		
-0.5597714,  0.3183714, 0.7916665,		
-0.2590142,  0.2301102, 0.7597134,		
-0.7580182,  0.5804536, 0.6936468,		
-0.2777158,  0.3438202, 0.6950511,		
0.4731839, -0.2248513, 0.7651013,		
0.5312460, -0.1991590, 0.7359806,		
-0.1663695,  0.3170776, 0.6480572,		
0.7478357, -0.2555715, 0.7105204),		
nrow=10, ncol=3)		
		
M2	<- matrix(	
c(1.0000000, 0.00000000,  0.000000000,  0.000000e+00,		
0.9071442, 0.04208416, -0.005683629,  0.000000e+00,		
0.8328467, 0.07604735, -0.012101182,  4.919680e-04,		
0.7643374, 0.09954517, -0.012344948, -6.355121e-05,		
0.7122795, 0.12114713, -0.016658244,  3.956041e-04,		
0.6621621, 0.13843375, -0.015044853, -5.779617e-04,		
0.6217507, 0.14719653, -0.013232027, -1.071607e-03,		
0.5818547, 0.15840180, -0.012182667, -1.541539e-03,		
0.5541449, 0.16511875, -0.010558783, -2.119333e-03,		
0.5209161, 0.17347345, -0.008961601, -2.680556e-03,		
0.5108935, 0.17432635, -0.009196419, -2.518288e-03),		
nrow=10, ncol=4)		

	bw<-beta.w
	mli<-mat_line
	daf_mean<-fail_data_mean
		
	if(bw<1.5)  {	
		a.0<- unname(M1a[mli,1]+ M1a[mli,2]*bw + M1a[mli,3]*bw^2)
	}else{	
		a.0<- unname(M1b[mli,1]+ M1b[mli,2]*bw + M1b[mli,3]*bw^2)
	}	
	b.0<-a.0/daf_mean *(M2[mli,1] + M2[mli,2]*log(a.0) + M2[mli,3]*log(a.0)^2 + M2[mli,4]*log(a.0)^2)	
	return(c(a.0,b.0))	
}		
