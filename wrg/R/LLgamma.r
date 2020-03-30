LLgamma<-function(x, s=NULL, alpha_g, beta_g)  {				
	suscomp<-0			
## Gbeta is also called the rate (1/scale), and is received as third arg to _gamma functions in R				
	failcomp<-sum(dgamma(x, alpha_g, beta_g, log=TRUE))			
	if(length(s)>0)  {			
		if(any(s<=0))  {		
			s2<-NULL	
			for(i in 1:length(s))  {	
				if(s[i]>0) {s2<-c(s2,s[i])}
			}	
			s<-s2	
			if(length(s)>0)  {	
				suscomp<-sum(pgamma(s, shape=alpha_g, rate=beta_g, log.p=TRUE))
			}	
		}else{		
			suscomp<-sum(pgamma(s, shape=alpha_g, rate=beta_g, log.p=TRUE))	
		}		
	}			
	value<-failcomp+suscomp			
	value			
}				
