# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:02:21 2019

@author: tengz
"""

##Input data
import numpy as np
import re
import sys
import numpy as np
import statsmodels.api as sm
import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
from scipy.stats import nbinom
from scipy.special import digamma
import statsmodels.sandbox as sms
from scipy.special import polygamma
###########
def fun_R_call(input_file, size_factor):

    def adj_loglikelihood(xVec,  lenSampleRna, X, y, mu, sign):    
        disp =  np.repeat(xVec, lenSampleRna)
        n = 1 / disp
        p = n / (n + mu)
        loglik = sum(nbinom.logpmf(y, n, p))
        diagVec = mu / (1 + np.dot(mu.transpose(), disp))
        diagWM = np.diagflat(diagVec)
        xtwx = np.dot(np.dot(np.transpose(X), diagWM), X)
        coxreid = 0.5 * np.log(np.linalg.det(xtwx))
        return (loglik - coxreid) * sign
    
    def calculate_varPrior(dispRaw, dispFitted, dispFittedIdx, varLogDispSamp):   
        logResidule = np.log(dispRaw[dispFittedIdx]) - np.log(dispFitted[dispFittedIdx])
        stdLogResidule = np.median(np.abs(logResidule - np.median(logResidule))) * 1.4826
        varLogResidule = stdLogResidule ** 2
        varPrior = varLogResidule - varLogDispSamp
        varPrior = max(varPrior, 0.1)
        return varPrior
    
	
    def calculate_logprior(disp, dispFitted, varPrior):	
        logprior = (np.log(disp) - np.log(dispFitted)) ** 2 / (2 * varPrior ** 2)
        return logprior
        
        

    def adj_loglikelihood_shrink(x, lenSampleRna, explanatory, response, yhat, dispFittedRna, varPriorRna, sign):	
        loglik_adj = adj_loglikelihood(x, lenSampleRna, explanatory, response, yhat, 1.0)
        logpriorRna = calculate_logprior(x, dispFittedRna, varPriorRna)
        loglik_adj_shrk = loglik_adj - logpriorRna
        return loglik_adj_shrk * sign
    
    
    
    
	
    fileNameCount = input_file   
    librarySizes = np.array(size_factor)
    with open(fileNameCount, 'r') as FileIn:
        header = np.array(FileIn.readline().strip().split(' '), dtype=str)
    geneIDs = np.loadtxt(fileNameCount, dtype=str, skiprows=1, usecols=(0,))
    geneIDs = geneIDs.reshape(geneIDs.size, 1)
    idxgene = header[1:int((len(header) + 1) / 2)]
    idxmethy = header[int((len(header) + 1) / 2):len(header)]
    idxGene = np.in1d(header, idxgene).nonzero()[0]
    idxMethy = np.in1d(header, idxmethy).nonzero()[0]
    countexpr = np.loadtxt(fileNameCount, dtype=int, skiprows=1, usecols=idxGene)
    countmethy = np.loadtxt(fileNameCount, dtype=float, skiprows=1, usecols=idxMethy)
    ## Caculate the raw disp

    num = len(geneIDs)

    lenSample = len(idxGene)
    muRaw = np.empty((num, lenSample))
    dispRaw = np.empty((num, 1))
    dispRaw.fill(np.nan)
    dispRawConv = dispRaw.copy()
    dispRawMthd = np.empty((num, 1), dtype='S10')
    dispRawMthd.fill('nan')
    dispInitial = 0.01
    sumCntCutoff = 10
    cntCutoff = sumCntCutoff
    responsea = countexpr
    explanatorya = countmethy

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print ('\r%i genes finished...' % i) ,
        if i+1 == num:
            print ('\r%i genes finished.' % num)

        if  sum(responsea[i, :] / librarySizes) >= cntCutoff:

            response =  responsea[i, :]
            explanatory =  np.hstack([np.ones(lenSample).reshape(lenSample,1), explanatorya[i,:].reshape(lenSample,1)])


            dispInitialRna  = dispInitial
            disp =  np.repeat(dispInitialRna, lenSample)

            mthd = 'SLSQP'
            j = 0
            while j < 10:
                try:
                    modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
                    result = modNB.fit()

                    dispBef = disp
                    x0 = dispBef[0]
                    yhat = result.mu
                    sign = -1.0

                    if mthd == 'SLSQP':
                        res = minimize(adj_loglikelihood, x0, args=(lenSample, explanatory, response, yhat, sign), method='SLSQP',  tol=1e-5)

                    if mthd == 'Nelder' or (mthd == 'SLSQP' and (not res.success or np.isnan(res.fun))):
                        if mthd == 'SLSQP':
                            j = 0
                            x0 =  dispInitialRna
                            dispBef =  np.repeat(dispInitialRna, lenSample)
                        res = minimize(adj_loglikelihood, x0, args=(lenSample, explanatory, response, yhat, sign), method='Nelder-Mead', tol=1e-5)
                        mthd = 'Nelder'

                    if res.success:
                        disp = np.repeat(res.x, lenSample)
                        if abs(np.log(disp[0]) - np.log(dispBef[0])) < 0.01:
                        
                            dispRaw[i]  = disp[0]
                            dispRawConv[i] = True
                            dispRawMthd[i] = mthd
                            muRaw[i, :] = yhat
                            break
                        elif j == 9:
                            
                            dispRaw[i] = disp[0]
                            dispRawConv[i] = False
                            dispRawMthd[i] = mthd
                         
                            muRaw[i, :] = yhat
                        else:
                            pass

                except sm.tools.sm_exceptions.PerfectSeparationError:
                    
                    dispRaw[i] = disp[0]
                    dispRawConv[i] = False
                    dispRawMthd[i] = mthd
                   
                    muRaw[i, :] = yhat

                j += 1

    if(len(np.where(np.isnan(dispRaw))[0])>0):
        na_num=len(np.where(np.isnan(dispRaw))[0])
        dispRaw[np.where(np.isnan(dispRaw))]= np.repeat(0.01,na_num)
    Idx = np.where(dispRaw<=2)[0]
    select_dispRaw = dispRaw[Idx]
    countMean = np.mean(countexpr / librarySizes,  axis=1)
    countMean = np.reshape(countMean, (countMean.size, 1))

    lowerBound = np.percentile(np.unique(select_dispRaw),  1)
    upperBound = np.percentile(np.unique(select_dispRaw), 99)
    idx = np.logical_and(select_dispRaw > lowerBound, select_dispRaw < upperBound).nonzero()[0]
    matrix = np.empty((idx.size, 2))
    matrix.fill(np.nan)
    matrix[:, 0] = 1 / countMean[idx].flatten()
    matrix[:, 1] = 1
    modGamma = sm.GLM(select_dispRaw[idx], matrix, family=sm.families.Gamma(sm.families.links.identity))
    result = modGamma.fit()
    Lambda = result.params
    dispFitted =select_dispRaw.copy()
    IDX = ~np.isnan(select_dispRaw)
    countMean = countMean[Idx]
    dispFitted[IDX] = Lambda[0] / countMean[IDX] + Lambda[1]
    ##Select Gene 
    select_countRNA = countexpr[Idx,]
    select_countmethy = countmethy[Idx,]
    gene_IDs = geneIDs[Idx]
    new_countRNA = select_countRNA[np.where(dispFitted>0)[0],]
    new_countmethy = select_countmethy[np.where(dispFitted>0)[0],]
    dispFitted = dispFitted[np.where(dispFitted>0)[0]]
    dispRaw = select_dispRaw[np.where(dispFitted>0)[0]]
    dispFittedIdx = ~np.isnan(dispRaw)
	
 

    num = len(dispRaw)
    muAdjRna = np.empty((num, idxGene.size))
    dispAdjRna = np.empty((num, 1))
    dispAdjRna.fill(np.nan)

    dispAdjConv = dispAdjRna.copy()
    dispAdjMthd = np.empty((num, 1), dtype='S10')
    dispAdjMthd.fill('nan')
    numSample = idxGene.size
    numCoef = 2
    varLogDispSamp = polygamma(1, (numSample - numCoef) / 2)
    dispFittedRna = dispFitted
    dispFittedRnaIdx = dispFittedIdx

    ## Caculate the prior disp
    varPriorRna = calculate_varPrior(dispRaw, dispFittedRna, dispFittedRnaIdx, varLogDispSamp)

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print ('\r%i genes finished...' % i) ,
        if i+1 == num:
            print ('\r%i genes finished.' % num)

        
        if not np.isnan(dispRaw[i]):
            response =  new_countRNA[i, :]
            explanatory = np.hstack([np.ones(lenSample).reshape(lenSample,1), new_countmethy[i,:].reshape(lenSample,1)])
            
            dispInitialRna  = dispInitial
            disp = np.repeat(dispInitialRna, lenSample)

            
            dispFitRna  = dispFittedRna[i]
            
            disprawRna  = dispRaw[i]
            mthd = 'SLSQP'
            j = 0
            while j < 10:
                try:
                    modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
                    result = modNB.fit()

                    dispBef = disp
                    x0 = dispBef[0]
                    yhat = result.mu
                    sign = -1.0

                    if mthd == 'SLSQP':
                        res = minimize(adj_loglikelihood_shrink, x0, args=(lenSample, explanatory, response, yhat,  dispFitRna, varPriorRna, sign), method='SLSQP', tol=1e-5)

                    if mthd == 'Nelder' or (mthd == 'SLSQP' and (not res.success or np.isnan(res.fun))):
                        if mthd == 'SLSQP':
                            j = 0
                            x0 =  dispInitialRna
                            dispBef =  np.repeat(dispInitialRna, lenSample)
                        res = minimize(adj_loglikelihood_shrink, x0, args=(lenSample, explanatory, response, yhat, dispFitRna, varPriorRna, sign), method='Nelder-Mead', tol=1e-5)
                        mthd = 'Nelder'

                    if res.success:
                        disp =  np.repeat(res.x, lenSample)
                        if abs(np.log(disp[0]) - np.log(dispBef[0])) < 0.01:
                            
                            dispAdjRna[i]  = disp[0]
                            dispAdjConv[i] = True
                            dispAdjMthd[i] = mthd
                        
                            muAdjRna[i, :] = yhat
                            break
                        elif j == 9:
                            
                            dispAdjRna[i]  = disp[0]
                            dispAdjConv[i] = False
                            dispAdjMthd[i] = mthd
                    
                            muAdjRna[i, :] = yhat
                        else:
                            pass

                except sm.tools.sm_exceptions.PerfectSeparationError:
                    
                    dispAdjRna[i] = disp[0]
                    dispAdjConv[i] = False
                    dispAdjMthd[i] = mthd
                    
                    muAdjRna[i, :] = yhat

                j += 1


    ## Calculate beta value
    if(len(np.where(np.isnan(dispAdjRna))[0])>0):
        na_num=len(np.where(np.isnan(dispAdjRna))[0])
        dispAdjRna[np.where(np.isnan(dispAdjRna))]= np.repeat(0.01,na_num)
    methy_betavalue = np.empty((num, 1))
    methy_betavalue.fill(np.nan)
    const_betavalue = np.empty((num, 1))
    const_betavalue.fill(np.nan)
    const_pvalue = np.empty((num, 1))
    const_pvalue.fill(np.nan)
    methy_pvalue = np.empty((num, 1))
    methy_pvalue.fill(np.nan)
    for i in range(num):
        sys.stdout.flush()

        if i % 50 == 0:
                print ('\r%i genes finished...' % i) ,
        if i+1 == num:
                print ('\r%i genes finished.' % num)
        if sum(new_countRNA[i, 0:int(lenSample/2)] / librarySizes[0:int(lenSample/2)]) >= 0 and sum(new_countRNA[i, int(lenSample/2):lenSample] / librarySizes[int(lenSample/2):lenSample]) >= 0:        
           if not np.isnan(dispAdjRna[i]):
                                           response =new_countRNA[i, :]
                                           if np.std(response)!=0:
                                               explanatory = np.hstack([np.ones(lenSample).reshape(lenSample,1), new_countmethy[i,:].reshape(lenSample,1)])
                                               disp = dispAdjRna[i]
                                               modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
                                               result0 = modNB.fit()
                                               const_betavalue[i] = result0.params[0]
                                               methy_betavalue[i] = result0.params[1]
                                               const_pvalue[i] = result0.pvalues[0]
                                               methy_pvalue[i] = result0.pvalues[1]
                                           
                                           if np.std(response)==0:
                                               i=i+1
                                               response=new_countRNA[i, :]
                                               explanatory = np.hstack([np.ones(lenSample).reshape(lenSample,1), new_countmethy[i,:].reshape(lenSample,1)])
                                               disp = dispAdjRna[i]
                                               modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
                                               result0 = modNB.fit()
                                               const_betavalue[i] = result0.params[0]
                                               methy_betavalue[i] = result0.params[1]
                                               const_pvalue[i] = result0.pvalues[0]
                                               methy_pvalue[i] = result0.pvalues[1]

    gene_ID = gene_IDs
    alpha = dispAdjRna
    const_betavalue = const_betavalue
    methy_betavalue = methy_betavalue
    const_pvalue = const_pvalue
    methy_pvalue = methy_pvalue

    Beta_Value = np.hstack([gene_ID, const_betavalue, methy_betavalue, const_pvalue, methy_pvalue, alpha])
    name = ['Gene_ID', 'Constant', 'Methylaiton', 'constant_p-value', 'methy_p-value', 'Alpha']
    out_result = pd.DataFrame(columns=name, data=Beta_Value) 
    return out_result
