#include <iostream>
#include "TH1.h"

Double_t bkgErrSqFromNorm(Double_t nEvtsBeforeNorm, Double_t normErrSq)
{
  return nEvtsBeforeNorm*nEvtsBeforeNorm*normErrSq;
}

Double_t bkgErrSqFromStats(Double_t norm, Double_t statErrBeforeNorm)
{
  return norm*norm*statErrBeforeNorm*statErrBeforeNorm;
}

Double_t bkgErrSqFromReweight(Double_t norm, Double_t reweightErrSq)
{
  return norm*norm*reweightErrSq;
}

Double_t normErrSqFromSearchStats(Double_t norm, Double_t nEvtsNormSearch)
{
  return norm*norm/nEvtsNormSearch;
}

Double_t normErrSqFromReweight(Double_t norm, Double_t nEvtsNormCtrlBeforeNorm, 
			       Double_t reweightErrSq)
{
  return norm*norm*reweightErrSq/(nEvtsNormCtrlBeforeNorm*nEvtsNormCtrlBeforeNorm);
}

Double_t normErrSqFromCtrlStats(Double_t norm, Double_t statErrSqBeforeNormCtrl, 
				Double_t nEvtsNormCtrlBeforeNorm)
{
  return norm*norm*statErrSqBeforeNormCtrl/(nEvtsNormCtrlBeforeNorm*nEvtsNormCtrlBeforeNorm);
}

Double_t normErrSq(const TH1* search, const TH1* beforeNorm, const TH1* reweightErr, 
		   Int_t normBinLow, Int_t normBinHigh, Double_t norm)
{
  Double_t retVal = -1.0;
  const Int_t nBinsSearch = search->GetNbinsX();
  const Int_t nBinsBeforeNorm = beforeNorm->GetNbinsX();
  const Int_t nBinsReweightErrSq = reweightErr->GetNbinsX();
  if (nBinsBeforeNorm != nBinsSearch) {
    cerr << "Error: beforeNorm->GetNbinsX() = " << nBinsBeforeNorm;
    cerr << " but search->GetNbinsX() = " << nBinsSearch << ".\n";
  }
  else if (nBinsBeforeNorm != nBinsReweightErrSq) {
    cerr << "Error: beforeNorm->GetNbinsX() = " << nBinsBeforeNorm;
    cerr << " but reweightErr->GetNbinsX() = " << nBinsReweightErrSq << ".\n";
  }
  else if (normBinLow > (nBinsSearch + 1)) {
    cerr << "Error: normBinLow = " << normBinLow << " but search->GetNbinsX() = " << nBinsSearch;
    cerr << ".\n";
  }
  else if (normBinHigh > (nBinsSearch + 1)) {
    cerr << "Error: normBinHigh = " << normBinHigh << " but search->GetNbinsX() = " << nBinsSearch;
    cerr << ".\n";
  }
  else {
    const Double_t nEvtsNormCtrlBeforeNorm = beforeNorm->Integral(normBinLow, normBinHigh);
    Double_t statErrSqBeforeNormCtrl = 0.0;
    Double_t reweightErrSqCtrl = 0.0;
    for (Int_t iBin = normBinLow; iBin <= normBinHigh; ++iBin) {
      const Double_t statErrSq = beforeNorm->GetBinError(iBin);
      const Double_t reweightErrSq = reweightErr->GetBinError(iBin);
      statErrSqBeforeNormCtrl+=statErrSq*statErrSq;
      reweightErrSqCtrl+=reweightErrSq*reweightErrSq;
    }
    retVal = normErrSqFromSearchStats(norm, search->Integral(normBinLow, normBinHigh)) + 
      normErrSqFromReweight(norm, nEvtsNormCtrlBeforeNorm, reweightErrSqCtrl) + 
      normErrSqFromCtrlStats(norm, statErrSqBeforeNormCtrl, nEvtsNormCtrlBeforeNorm);
  }
  return retVal;
}

Double_t bkgErrSq(const TH1* beforeNorm, const TH1* reweightErr, Int_t bin, Double_t norm, 
		  Double_t normErrSq)
{
  Double_t retVal = -1.0;
  const Int_t nBinsBeforeNorm = beforeNorm->GetNbinsX();
  const Int_t nBinsReweightErr = reweightErr->GetNbinsX();
  if (nBinsBeforeNorm != nBinsReweightErr) {
    cerr << "Error: beforeNorm->GetNbinsX() = " << nBinsBeforeNorm;
    cerr << " but reweightErr->GetNbinsX() = " << nBinsReweightErr << ".\n";
  }
  else if (bin > (nBinsBeforeNorm + 1)) {
    cerr << "Error: bin = " << bin << " but beforeNorm->GetNbinsX() = " << nBinsBeforeNorm;
    cerr << ".\n";
  }
  else {
    const Double_t binReweightErr = reweightErr->GetBinError(bin);
    retVal = bkgErrSqFromNorm(beforeNorm->GetBinContent(bin), normErrSq) + 
      bkgErrSqFromStats(norm, beforeNorm->GetBinError(bin)) + 
      bkgErrSqFromReweight(norm, binReweightErr*binReweightErr);
  }
  return retVal;
}
