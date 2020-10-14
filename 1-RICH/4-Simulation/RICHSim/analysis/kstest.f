SUBROUTINE ks2d1s(x1,y1,n1,quadvl,d1,prob) 
    INTEGER n1
    REAL d1,prob,x1(n1),y1(n1)
    EXTERNAL quadvl
C USES pearsn,probks,quadct,quadvl Two-dimensional Kolmogorov-Smirnov test of one sample against a model. Given the x and y coordinates of n1 data points in arrays x1(1:n1) and y1(1:n1), and given a user-supplied function quadvl that exemplifies the model, this routine returns the two- dimensional K-S statistic as d1, and its significance level as prob. Small values of prob show that the sample is significantly different from the model. Note that the test is slightly distribution-dependent, so prob is only an estimate.
    INTEGER j
    REAL dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,rr,sqen,probks 
    d1=0.0
    do j=1,n1 
C Loop over the data points. 
        call quadct(x1(j),y1(j),x1,y1,n1,fa,fb,fc,fd)
        call quadvl(x1(j),y1(j),ga,gb,gc,gd) 
        d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
C For both the sample and the model, the distribution is integrated in each of four quad-rants, and the maximum difference is saved.
    enddo
    call pearsn(x1,y1,n1,r1,dum,dumm)
C Get the linear correlation coefficient r1. 
    sqen=sqrt(float(n1))
    rr=sqrt(1.0-r1**2)
C Estimate the probability using the K-S probability function probks. 
    prob=probks(d1*sqen/(1.0+rr*(0.25-0.75/sqen)))
    return
END

SUBROUTINE quadct(x,y,xx,yy,nn,fa,fb,fc,fd) INTEGER nn
    REAL fa,fb,fc,fd,x,y,xx(nn),yy(nn)
C Given an origin (x, y), and an array of nn points with coordinates xx and yy, count how many of them are in each quadrant around the origin, and return the normalized frac- tions. Quadrants are labeled alphabetically, counterclockwise from the upper right. Used by ks2d1s and ks2d2s.
    INTEGER k,na,nb,nc,nd 
    REAL ff
    na=0
    nb=0
    nc=0 
    nd=0 
    do 11 k=1,nn
        if(yy(k).gt.y)then
            if(xx(k).gt.x)then 
                na=na+1
            else 
                nb=nb+1
            endif
        else
            if(xx(k).gt.x)then 
                nd=nd+1
            else 
                nc=nc+1
            endif 
        endif
    enddo 11 
    ff=1.0/nn 
    fa=ff*na 
    fb=ff*nb
    fc=ff*nc
    fd=ff*nd
    return
END

SUBROUTINE quadvl(x,y,fa,fb,fc,fd) 
    REAL fa,fb,fc,fd,x,y
C This is a sample of a user-supplied routine to be used with ks2d1s. In this case, the model distribution is uniform inside the square −1 < x < 1, −1 < y < 1. In general this routine should return, for any point (x,y), the fraction of the total distribution in each of the four quadrants around that point. The fractions, fa, fb, fc, and fd, must add up to 1. Quadrants are alphabetical, counterclockwise from the upper right.
    REAL qa,qb,qc,qd 
    qa=min(2.,max(0.,1.-x)) 
    qb=min(2.,max(0.,1.-y)) 
    qc=min(2.,max(0.,x+1.)) 
    qd=min(2.,max(0.,y+1.)) 
    fa=0.25*qa*qb 
    fb=0.25*qb*qc 
    fc=0.25*qc*qd 
    fd=0.25*qd*qa
    return 
END

SUBROUTINE ks2d2s(x1,y1,n1,x2,y2,n2,d,prob) 
    INTEGER n1,n2
    REAL d,prob,x1(n1),x2(n2),y1(n1),y2(n2)
C USESpearsn,probks,quadct Two-dimensional Kolmogorov-Smirnov test on two samples. Given the x and y coordinates of the first sample as n1 values in arrays x1(1:n1) and y1(1:n1), and likewise for the second sample, n2 values in arrays x2 and y2, this routine returns the two-dimensional, two- sample K-S statistic as d, and its significance level as prob. Small values of prob show that the two samples are significantly different. Note that the test is slightly distribution- dependent, so prob is only an estimate.
    INTEGER j
    REAL d1,d2,dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,r2,rr,sqen,probks
    d1=0.0
    do 11 j=1,n1 
CFirst, use points in the first sample as origins. 
        call quadct(x1(j),y1(j),x1,y1,n1,fa,fb,fc,fd)
        call quadct(x1(j),y1(j),x2,y2,n2,ga,gb,gc,gd) 
        d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
    enddo 11
    d2=0.0
    do 12 j=1,n2 
C Then, use points in the second sample as origins.
        call quadct(x2(j),y2(j),x1,y1,n1,fa,fb,fc,fd)
        call quadct(x2(j),y2(j),x2,y2,n2,ga,gb,gc,gd) 
        d2=max(d2,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
    enddo 12

    d=0.5*(d1+d2) 
C Average the K-S statistics. 
    sqen=sqrt(float(n1)*float(n2)/float(n1+n2))
    call pearsn(x1,y1,n1,r1,dum,dumm) 
    call pearsn(x2,y2,n2,r2,dum,dumm)
C Get the linear correlation coefficient for each sample. 
    rr=sqrt(1.0-0.5*(r1**2+r2**2))
C Estimate the probability using the K-S probability function probks. 
    prob=probks(d*sqen/(1.0+rr*(0.25-0.75/sqen)))
    return
END