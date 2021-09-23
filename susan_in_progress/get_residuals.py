def linregress_residuals(xdata,ydata):
    """
    This function performs a linear regression and then gets the residuals
    
    Args:
        xdata (array-like): The x data
        ydata (array-like): The y data
        
    Returns:
        residuals: the residuals of the regression
        slope: the slope of regression line
        intercept: intercept of the regression line
        rvalue: correlation coeffficient
        pvalue: two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero.
        stderr: standard error of the estimated gradient
    
    Author: SMM
    """
    
    from scipy import stats
    
    # Get the regression
    (m,b,r,pvalue,stderr)=stats.linregress(xdata,ydata)   
    #print("In regress1, m is: "+str(m))    
    # get the residuals
    residuals = np.copy(xdata)
    for idx,x in enumerate(xdata):
        yfit = m*x+b
        residuals[idx] = yfit-ydata[idx]
        
    return (residuals,m,b,r,pvalue,stderr) 

